use std::collections::HashMap;
use std::path::PathBuf;
use std::str;
use std::sync::Arc;

use anyhow::{Context, Result};
use bio::stats::{bayesian, LogProb};
use crossbeam::channel::{Receiver, Sender};
use derive_builder::Builder;
use itertools::Itertools;
use rust_htslib::bcf::{self, Read};

use crate::calling::variants::preprocessing::{
    read_observations, remove_observation_header_entries, OBSERVATION_FORMAT_VERSION,
};
use crate::calling::variants::SampleInfo;
use crate::calling::variants::{
    chrom, event_tag_name, Call, CallBuilder, SampleInfoBuilder, VariantBuilder,
};
use crate::errors;
use crate::grammar;
use crate::utils;
use crate::utils::worker_pool;
use crate::utils::worker_pool::Orderable;
use crate::variants::evidence::observation::Observation;
use crate::variants::model;
use crate::variants::model::modes::generic::{
    self, GenericLikelihood, GenericModelBuilder, GenericPosterior,
};
use crate::variants::model::Contamination;
use crate::variants::model::{AlleleFreq, StrandBias};

pub(crate) type AlleleFreqCombination = Vec<model::likelihood::Event>;

pub(crate) type Model<Pr> =
    bayesian::Model<GenericLikelihood, Pr, GenericPosterior, generic::Cache>;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct Caller<Pr>
where
    Pr: bayesian::model::Prior,
{
    samplenames: grammar::SampleInfo<String>,
    observations: grammar::SampleInfo<PathBuf>,
    scenario: grammar::Scenario,
    outbcf: Option<PathBuf>,
    contaminations: grammar::SampleInfo<Option<Contamination>>,
    resolutions: grammar::SampleInfo<usize>,
    prior: Pr,
    threads: usize,
    buffer_capacity: usize,
}

impl<Pr> Caller<Pr>
where
    Pr: bayesian::model::Prior<Event = AlleleFreqCombination>
        + model::modes::UniverseDrivenPrior
        + Clone
        + Default
        + Send
        + Sync,
{
    pub(crate) fn n_samples(&self) -> usize {
        self.samplenames.len()
    }

    pub(crate) fn header(&self) -> Result<bcf::Header> {
        let mut header = bcf::Header::from_template(
            bcf::Reader::from_path(self.observations.first().as_ref().unwrap())?.header(),
        );

        remove_observation_header_entries(&mut header);

        // register samples
        for sample_name in self.samplenames.iter() {
            header.push_sample(sample_name.as_bytes());
        }

        // register events
        for event in self.scenario.events().keys() {
            header.push_record(
                format!(
                    "##INFO=<ID={},Number=A,Type=Float,\
                     Description=\"Posterior probability for event {} (PHRED)\">",
                    event_tag_name(event),
                    event
                )
                .as_bytes(),
            );
        }
        header.push_record(
            b"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,\
             Description=\"Posterior probability for strand bias artifact (PHRED)\">",
        );
        header.push_record(
            b"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,\
             Description=\"Posterior probability for not having a variant (PHRED)\">",
        );

        // register sample specific tags
        header.push_record(
            b"##FORMAT=<ID=DP,Number=A,Type=Integer,\
              Description=\"Expected sequencing depth, while considering mapping uncertainty\">",
        );
        header.push_record(
            b"##FORMAT=<ID=OBS,Number=A,Type=String,\
              Description=\"Posterior odds for alt allele of each fragment as Kass Raftery \
              scores: N=none, B=barely, P=positive, S=strong, V=very strong (lower case if \
              probability for correct mapping of fragment is <95%)\">",
        );
        header.push_record(
            b"##FORMAT=<ID=AF,Number=A,Type=Float,\
              Description=\"Maximum a posteriori probability estimate of allele frequency\">",
        );
        header.push_record(
            b"##FORMAT=<ID=SB,Number=A,Type=String,\
              Description=\"Strand bias estimate: + indicates that ALT allele is associated with \
              forward strand, - indicates that ALT allele is associated with reverse strand, \
              - indicates no strand bias.\">",
        );

        Ok(header)
    }

    pub(crate) fn writer(&self) -> Result<bcf::Writer> {
        let header = self.header();

        Ok(if let Some(ref path) = self.outbcf {
            bcf::Writer::from_path(path, &header.as_ref().unwrap(), false, bcf::Format::BCF)
                .context(format!("Unable to write BCF to {}.", path.display()))?
        } else {
            bcf::Writer::from_stdout(&header.as_ref().unwrap(), false, bcf::Format::BCF)
                .context("Unable to write BCF to STDOUT.")?
        })
    }

    fn model(&self) -> Model<Pr> {
        GenericModelBuilder::default()
            // TODO allow to define prior in the grammar
            .prior(self.prior.clone())
            .contaminations(self.contaminations.clone())
            .resolutions(self.resolutions.clone())
            .build()
            .unwrap()
    }

    fn observations(&self) -> Result<grammar::SampleInfo<bcf::Reader>> {
        let mut observations = grammar::SampleInfo::default();
        for path in self.observations.iter() {
            observations.push(bcf::Reader::from_path(path)?);
        }
        Ok(observations)
    }

    pub(crate) fn call(&self) -> Result<()> {
        // Configure worker pool:
        let preprocessor = |sender: Sender<Vec<WorkItem>>,
                            buffer_guard: Arc<utils::worker_pool::BufferGuard>|
         -> Result<()> {
            let mut observations = self.observations()?;

            // Check observation format.
            for obs_reader in observations.iter() {
                let mut valid = false;
                for record in obs_reader.header().header_records() {
                    if let bcf::HeaderRecord::Generic { key, value } = record {
                        if key == "varlociraptor_observation_format_version"
                            && value == OBSERVATION_FORMAT_VERSION
                        {
                            valid = true;
                        }
                    }
                }
                if !valid {
                    return Err(errors::Error::InvalidObservationFormat.into());
                }
            }

            let mut work_items = Vec::new();
            let mut last_bnd_event = None;
            let mut i = 0;
            loop {
                let mut records = observations.map(|reader| reader.empty_record());
                let mut eof = Vec::new();
                for (reader, record) in observations.iter_mut().zip(records.iter_mut()) {
                    eof.push(!reader.read(record)?);
                }

                if eof.iter().all(|v| *v) {
                    return Ok(());
                } else if !eof.iter().all(|v| !v) {
                    // only some are EOF, this is an error
                    return Err(errors::Error::InconsistentObservations.into());
                }

                // ensure that all observation BCFs contain exactly the same calls
                let first_record = records.first().unwrap();
                let current_rid = first_record.rid();
                let current_pos = first_record.pos();
                let current_alleles = first_record.alleles();
                for record in &records[1..] {
                    if record.rid() != current_rid
                        || record.pos() != current_pos
                        || record.alleles() != current_alleles
                    {
                        return Err(errors::Error::InconsistentObservations.into());
                    }
                }

                let work_item =
                    self.preprocess_record(&mut records, i, &observations, &last_bnd_event)?;
                let do_call = match (&work_item.bnd_event, &last_bnd_event) {
                    (Some(current), Some(last)) if current == last => {
                        // Same breakend event, do not call yet.
                        false
                    }
                    (Some(current), _) => {
                        // Other or no breakend event in current record, hence call.
                        last_bnd_event = Some(current.to_owned());
                        true
                    }
                    (None, _) => {
                        // No breakend event at all, hence call.
                        last_bnd_event = None;
                        true
                    }
                };
                work_items.push(work_item);

                if do_call {
                    // send vector of call infos
                    let tosend = work_items;
                    sender.send(tosend).unwrap();
                    // clear vector
                    work_items = Vec::new();
                }
                buffer_guard.wait_for_free();

                i += 1;
            }
        };

        let mut workers = Vec::new();
        for _ in 0..self.threads {
            workers.push(
                |receiver: Receiver<Vec<WorkItem>>, sender: Sender<Box<WorkItem>>| -> Result<()> {
                    let mut model = self.model();
                    let mut events = Vec::new();
                    let mut last_rid = None;
                    let mut breakend_result = None;
                    for work_items in receiver {
                        for mut work_item in work_items {
                            let contig = str::from_utf8(work_item.call.chrom()).unwrap();
                            self.configure_model(
                                work_item.rid,
                                last_rid,
                                &mut model,
                                &mut events,
                                contig,
                            )?;
                            last_rid = Some(work_item.rid);

                            self.call_record(&mut work_item, &model, &events, &mut breakend_result);
                            sender.send(Box::new(work_item)).unwrap();
                        }
                    }
                    Ok(())
                },
            );
        }

        let postprocessor = |receiver: Receiver<Box<WorkItem>>| -> Result<()> {
            let mut bcf_writer = self.writer()?;
            for work_item in receiver {
                work_item.call.write_final_record(&mut bcf_writer)?;

                if work_item.index % 100 == 0 {
                    info!("{} records processed.", work_item.index);
                }
            }

            Ok(())
        };

        worker_pool(
            preprocessor,
            workers.iter(),
            postprocessor,
            self.buffer_capacity,
        )
    }

    fn preprocess_record(
        &self,
        records: &mut grammar::SampleInfo<bcf::Record>,
        index: usize,
        observations: &grammar::SampleInfo<bcf::Reader>,
        last_breakend_event: &Option<Vec<u8>>,
    ) -> Result<WorkItem> {
        let (call, snv, bnd_event, rid) = {
            let first_record = records
                .first_mut()
                .expect("bug: there must be at least one record");
            let start = first_record.pos() as u64;
            let chrom = chrom(
                &observations
                    .first()
                    .expect("bug: there must be at least one observation reader"),
                first_record,
            );

            let call = CallBuilder::default()
                .chrom(chrom.to_owned())
                .pos(start)
                .id({
                    let id = first_record.id();
                    if id == b"." {
                        None
                    } else {
                        Some(id)
                    }
                })
                .variants(Vec::new())
                .record(first_record)?
                .build()
                .unwrap();

            // store information about SNV for special handling in posterior computation (variant selection operations)
            let snv = {
                let alleles = first_record.alleles();
                if alleles[0].len() == 1 && alleles[1].len() == 1 {
                    Some(model::modes::generic::SNV::new(
                        alleles[0][0],
                        alleles[1][0],
                    ))
                } else {
                    None
                }
            };

            let bnd_event = if utils::is_bnd(first_record)? {
                Some(utils::info_tag_event(first_record)?.unwrap().to_owned())
            } else {
                None
            };

            let rid = first_record
                .rid()
                .ok_or_else(|| errors::Error::RecordMissingChrom { i: index + 1 })?;

            (call, snv, bnd_event, rid)
        };

        let mut variant_builder = VariantBuilder::default();
        variant_builder.record(records.first_mut().unwrap())?;

        let mut work_item = WorkItem {
            rid,
            call,
            pileups: None,
            snv,
            bnd_event,
            variant_builder,
            index,
        };

        if let Some(ref event) = work_item.bnd_event {
            if let Some(ref last_event) = last_breakend_event {
                if last_event == event {
                    // METHOD: Another breakend in the same event was already processed, hence, we will just copy the
                    // results (no pileup needed). This works because preprocessing ensures that all breakends of one event appear
                    // "en block".
                    return Ok(work_item);
                }
            }
        }

        let mut pileups = Vec::new();
        for record in records.iter_mut() {
            pileups.push(read_observations(record)?);
        }

        // obtain pileups
        work_item.pileups = Some(pileups);

        Ok(work_item)
    }

    fn configure_model(
        &self,
        current_rid: u32,
        rid: Option<u32>,
        model: &mut Model<Pr>,
        events: &mut Vec<model::Event>,
        contig: &str,
    ) -> Result<()> {
        if !rid.map_or(false, |rid: u32| current_rid == rid) {
            // rid is not the same as before, obtain event universe
            // clear old events
            events.clear();

            // register absent event
            events.push(model::Event {
                name: "absent".to_owned(),
                vafs: grammar::VAFTree::absent(self.n_samples()),
                strand_bias: StrandBias::None,
            });

            // add events from scenario
            for (event_name, vaftree) in self.scenario.vaftrees(contig)? {
                events.push(model::Event {
                    name: event_name.clone(),
                    vafs: vaftree.clone(),
                    strand_bias: StrandBias::None,
                });
                events.push(model::Event {
                    name: event_name.clone(),
                    vafs: vaftree.clone(),
                    strand_bias: StrandBias::Forward,
                });
                events.push(model::Event {
                    name: event_name,
                    vafs: vaftree,
                    strand_bias: StrandBias::Reverse,
                });
            }

            // update prior to the VAF universe of the current chromosome
            let mut vaf_universes = self.scenario.sample_info();
            for (sample_name, sample) in self.scenario.samples().iter() {
                let universe = sample.contig_universe(&contig)?;
                vaf_universes = vaf_universes.push(sample_name, universe.to_owned());
            }

            model.prior_mut().set_universe(vaf_universes.build());
        }

        Ok(())
    }

    fn call_record(
        &self,
        work_item: &mut WorkItem,
        model: &Model<Pr>,
        event_universe: &[model::Event],
        breakend_result: &mut Option<BreakendResult>,
    ) {
        // In the else case, the variant builder can be filled with info from a previous breakend.
        if work_item.pileups.is_some() {
            let data = model::modes::generic::Data::new(
                work_item.pileups.take().unwrap(),
                work_item.snv.clone(),
            );

            // Compute probabilities for given events.
            let m = model.compute(event_universe.iter().cloned(), &data);

            // add calling results
            let mut event_probs: HashMap<String, LogProb> = event_universe
                .iter()
                .filter_map(|event| {
                    if event.is_artifact() {
                        None
                    } else {
                        Some((event.name.clone(), m.posterior(event).unwrap()))
                    }
                })
                .collect();
            // generate artifact event
            event_probs.insert(
                "artifact".to_owned(),
                LogProb::ln_sum_exp(
                    &event_universe
                        .iter()
                        .filter_map(|event| {
                            if event.is_artifact() {
                                Some(m.posterior(event).unwrap())
                            } else {
                                None
                            }
                        })
                        .collect_vec(),
                ),
            );
            work_item.variant_builder.event_probs(Some(event_probs));

            // add sample specific information
            work_item.variant_builder.sample_info(
                if let Some(map_estimates) = m.maximum_posterior() {
                    Some(
                        data.into_pileups()
                            .into_iter()
                            .zip(map_estimates.iter())
                            .map(|(pileup, estimate)| {
                                let mut sample_builder = SampleInfoBuilder::default();
                                sample_builder.observations(pileup);
                                match estimate {
                                    model::likelihood::Event { strand_bias, .. }
                                        if strand_bias.is_some() =>
                                    {
                                        sample_builder
                                            .allelefreq_estimate(AlleleFreq(0.0))
                                            .strand_bias(*strand_bias);
                                    }
                                    model::likelihood::Event { allele_freq, .. } => {
                                        sample_builder
                                            .allelefreq_estimate(*allele_freq)
                                            .strand_bias(StrandBias::None);
                                    }
                                };
                                sample_builder.build().unwrap()
                            })
                            .collect_vec(),
                    )
                } else {
                    // no observations
                    Some(Vec::new())
                },
            );
        } else if let Some(result) = breakend_result {
            // Take sample info and event probs from previous breakend.
            work_item
                .variant_builder
                .event_probs(Some(result.event_probs.clone()));
            work_item
                .variant_builder
                .sample_info(Some(result.sample_info.clone()));
        } else {
            unreachable!();
        }

        let variant = work_item.variant_builder.build().unwrap();

        if let Some(ref event) = work_item.bnd_event {
            if breakend_result
                .as_ref()
                .map_or(false, |result| &result.event != event)
            {
                // Store breakend results for next breakend of the same event.
                breakend_result.replace(BreakendResult {
                    event: event.to_owned(),
                    event_probs: variant.event_probs().as_ref().unwrap().clone(),
                    sample_info: variant.sample_info().as_ref().unwrap().clone(),
                });
            }
        } else {
            breakend_result.take();
        }

        work_item.call.variants.push(variant);
    }
}

struct BreakendResult {
    event: Vec<u8>,
    event_probs: HashMap<String, LogProb>,
    sample_info: Vec<SampleInfo>,
}

struct WorkItem {
    rid: u32,
    call: Call,
    variant_builder: VariantBuilder,
    pileups: Option<Vec<Vec<Observation>>>,
    snv: Option<model::modes::generic::SNV>,
    bnd_event: Option<Vec<u8>>,
    index: usize,
}

impl Orderable for WorkItem {
    fn index(&self) -> usize {
        self.index
    }
}

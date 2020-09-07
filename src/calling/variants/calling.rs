use std::collections::HashMap;
use std::path::PathBuf;
use std::str;
use std::sync::{Mutex, RwLock};

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
use crate::variants::evidence::observation::Observation;
use crate::variants::model;
use crate::variants::model::modes::generic::{
    self, GenericLikelihood, GenericModelBuilder, GenericPosterior,
};
use crate::variants::model::Contamination;
use crate::variants::model::{
    bias::read_orientation_bias::ReadOrientationBias, bias::strand_bias::StrandBias, bias::Biases,
    bias::BiasesBuilder, AlleleFreq,
};
use crate::variants::types::breakends::BreakendIndex;

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
    breakend_index: BreakendIndex,
    #[builder(default)]
    breakend_results: RwLock<HashMap<Vec<u8>, BreakendResult>>,
    threads: usize,
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
              Description=\"Summary of observations. Each entry is encoded as CBTSO, with C being a count, \
              B being the posterior odds for the alt allele (see below), T being the type of alignment, encoded \
              as s=single end and p=paired end, S being the strand that supports the observation (+, -, or * for both), \
              and O being the read orientation (> = F1R2, < = F2R1, * = unknown, ! = non standard, e.g. R1F2). \
              Posterior odds for alt allele of each fragment are given as Kass Raftery \
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
              - indicates no strand bias. Strand bias is indicative for systematic sequencing \
              errors. Probability for strand bias is captured by the ARTIFACT event (PROB_ARTIFACT).\">",
        );
        header.push_record(
            b"##FORMAT=<ID=ROB,Number=A,Type=String,\
              Description=\"Read orientation bias estimate: > indicates that ALT allele is associated with \
              F1R2 orientation, < indicates that ALT allele is associated with F2R1 orientation, \
              - indicates no read orientation bias. Read orientation bias is indicative of Guanin \
              oxidation artifacts. Probability for read orientation bias is captured by the ARTIFACT \
              event (PROB_ARTIFACT).\">",
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
        let preprocessor = |sender: Sender<WorkItem>| -> Result<()> {
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

                let work_item = self.preprocess_record(&mut records, i, &observations)?;

                sender.send(work_item).unwrap();

                i += 1;
            }
        };

        let mut workers = Vec::new();
        for _ in 0..self.threads {
            workers.push(
                |receiver: Receiver<WorkItem>, sender: Sender<WorkItem>| -> Result<()> {
                    let mut model = self.model();
                    // For SNVs and MNVs we need a special model as here read orientation bias needs to be considered.
                    let mut read_orientation_bias_model = self.model();
                    let mut events = Vec::new();
                    let mut last_rid = None;
                    let mut last_read_orientation_bias_rid = None;
                    for mut work_item in receiver {
                        let contig = str::from_utf8(work_item.call.chrom()).unwrap();
                        let _model;
                        let _last_rid;

                        if work_item.check_read_orientation_bias {
                            _model = &mut read_orientation_bias_model;
                            _last_rid = last_read_orientation_bias_rid;
                            last_read_orientation_bias_rid = Some(work_item.rid);
                        } else {
                            _model = &mut model;
                            _last_rid = last_rid;
                            last_rid = Some(work_item.rid);
                        }

                        self.configure_model(
                            work_item.rid,
                            _last_rid,
                            _model,
                            &mut events,
                            contig,
                            work_item.check_read_orientation_bias,
                        )?;

                        self.call_record(&mut work_item, _model, &events);
                        sender.send(work_item).unwrap();
                    }
                    Ok(())
                },
            );
        }

        let postprocessor = |receiver: Receiver<WorkItem>| -> Result<()> {
            let mut bcf_writer = self.writer()?;
            let mut processed = 0;
            for work_item in receiver {
                work_item.call.write_final_record(&mut bcf_writer)?;
                processed += 1;

                if processed % 100 == 0 {
                    info!("{} records processed.", processed);
                }
            }

            Ok(())
        };

        worker_pool(preprocessor, workers.iter(), postprocessor)
    }

    fn preprocess_record(
        &self,
        records: &mut grammar::SampleInfo<bcf::Record>,
        index: usize,
        observations: &grammar::SampleInfo<bcf::Reader>,
    ) -> Result<WorkItem> {
        let (call, snv, bnd_event, rid, is_snv_or_mnv) = {
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
            let snv;
            let is_snv_or_mnv;
            {
                let alleles = first_record.alleles();
                if alleles[0].len() == 1 && alleles[1].len() == 1 {
                    is_snv_or_mnv = true;
                    snv = Some(model::modes::generic::SNV::new(
                        alleles[0][0],
                        alleles[1][0],
                    ));
                } else if alleles[0].len() == alleles[1].len() {
                    is_snv_or_mnv = true;
                    snv = None;
                } else {
                    is_snv_or_mnv = false;
                    snv = None;
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

            (call, snv, bnd_event, rid, is_snv_or_mnv)
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
            check_read_orientation_bias: is_snv_or_mnv,
        };

        if let Some(ref event) = work_item.bnd_event {
            if self.breakend_results.read().unwrap().contains_key(event) {
                // METHOD: Another breakend in the same event was already processed, hence, we will just copy the
                // results (no pileup needed).
                return Ok(work_item);
            }
        }

        // obtain pileups
        let mut paired_end = false;
        let mut pileups = Vec::new();
        for record in records.iter_mut() {
            let pileup = read_observations(record)?;
            paired_end |= pileup.iter().any(|obs| obs.is_paired());
            pileups.push(pileup);
        }

        work_item.pileups = Some(pileups);
        // Only check for read orientation bias if there is at least one paired end read.
        work_item.check_read_orientation_bias &= paired_end;

        Ok(work_item)
    }

    fn configure_model(
        &self,
        current_rid: u32,
        rid: Option<u32>,
        model: &mut Model<Pr>,
        events: &mut Vec<model::Event>,
        contig: &str,
        consider_read_orientation_bias: bool,
    ) -> Result<()> {
        if !rid.map_or(false, |rid: u32| current_rid == rid) {
            // rid is not the same as before, obtain event universe
            // clear old events
            events.clear();

            // register absent event
            events.push(model::Event {
                name: "absent".to_owned(),
                vafs: grammar::VAFTree::absent(self.n_samples()),
                biases: Biases::none(),
            });

            // add events from scenario
            for (event_name, vaftree) in self.scenario.vaftrees(contig)? {
                events.push(model::Event {
                    name: event_name.clone(),
                    vafs: vaftree.clone(),
                    biases: Biases::none(),
                });
                // Corresponding biased events.
                for biases in Biases::all_artifact_combinations(consider_read_orientation_bias) {
                    events.push(model::Event {
                        name: event_name.clone(),
                        vafs: vaftree.clone(),
                        biases,
                    });
                }
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
    ) {
        if let Some(ref bnd_event) = work_item.bnd_event {
            if let Some(result) = self.breakend_results.read().unwrap().get(bnd_event) {
                // Take sample info and event probs from previous breakend.
                work_item
                    .variant_builder
                    .event_probs(Some(result.event_probs.clone()));
                work_item
                    .variant_builder
                    .sample_info(Some(result.sample_info.clone()));

                let variant = work_item.variant_builder.build().unwrap();
                work_item.call.variants.push(variant);

                return;
            }
        }

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
                                    model::likelihood::Event { biases, .. }
                                        if biases.is_artifact() =>
                                    {
                                        sample_builder
                                            .allelefreq_estimate(AlleleFreq(0.0))
                                            .biases(biases.clone());
                                    }
                                    model::likelihood::Event { allele_freq, .. } => {
                                        sample_builder
                                            .allelefreq_estimate(*allele_freq)
                                            .biases(Biases::none());
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
        } else {
            unreachable!();
        }

        let variant = work_item.variant_builder.build().unwrap();

        if let Some(ref event) = work_item.bnd_event {
            if self.breakend_index.last_record_index(event).unwrap() == work_item.index {
                // METHOD: last index, hence clear result
                self.breakend_results.write().unwrap().remove(event);
            } else {
                // METHOD: store breakend group result for next breakend of this group
                self.breakend_results.write().unwrap().insert(
                    event.to_owned(),
                    BreakendResult {
                        event: event.to_owned(),
                        event_probs: variant.event_probs().as_ref().unwrap().clone(),
                        sample_info: variant.sample_info().as_ref().unwrap().clone(),
                    },
                );
            }
        }

        work_item.call.variants.push(variant);
    }
}

#[derive(Default)]
pub(crate) struct BreakendResult {
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
    check_read_orientation_bias: bool,
}

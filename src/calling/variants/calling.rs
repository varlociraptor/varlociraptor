use std::collections::HashMap;
use std::error::Error;
use std::path::Path;
use std::str;

use bio::stats::{bayesian, LogProb};
use derive_builder::Builder;
use itertools::Itertools;
use rust_htslib::bcf::{self, Read};

use crate::calling::variants::preprocessing::{
    read_observations, remove_observation_header_entries, OBSERVATION_FORMAT_VERSION,
};
use crate::calling::variants::{
    chrom, event_tag_name, Call, CallBuilder, SampleInfoBuilder, VariantBuilder,
};
use crate::errors;
use crate::grammar;
use crate::model;
use crate::model::{AlleleFreq, StrandBias};

pub type AlleleFreqCombination = Vec<model::likelihood::Event>;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller<L, Pr, Po, ModelPayload>
where
    L: bayesian::model::Likelihood<ModelPayload>,
    Pr: bayesian::model::Prior,
    Po: bayesian::model::Posterior<Event = model::Event>,
    ModelPayload: Default,
{
    samplenames: grammar::SampleInfo<String>,
    observations: grammar::SampleInfo<bcf::Reader>,
    scenario: grammar::Scenario,
    #[builder(private)]
    bcf_writer: bcf::Writer,
    model: bayesian::Model<L, Pr, Po, ModelPayload>,
}

impl<L, Pr, Po, ModelPayload> CallerBuilder<L, Pr, Po, ModelPayload>
where
    L: bayesian::model::Likelihood<ModelPayload>,
    Pr: bayesian::model::Prior,
    Po: bayesian::model::Posterior<Event = model::Event>,
    ModelPayload: Default,
{
    pub fn outbcf<P: AsRef<Path>>(self, path: Option<P>) -> Result<Self, Box<dyn Error>> {
        let mut header = bcf::Header::from_template(
            self.observations
                .as_ref()
                .expect("bug: observations() must be called bfore outbcf()")
                .first()
                .unwrap()
                .header(),
        );

        remove_observation_header_entries(&mut header);

        // register samples
        for sample_name in self
            .samplenames
            .as_ref()
            .expect(".samples() has to be called before .outbcf()")
            .iter()
        {
            header.push_sample(sample_name.as_bytes());
        }

        // register events
        for event in self
            .scenario
            .as_ref()
            .expect(".scenario() has to be called before .outbcf()")
            .events()
            .keys()
        {
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

        let writer = if let Some(path) = path {
            bcf::Writer::from_path(path, &header, false, bcf::Format::BCF)?
        } else {
            bcf::Writer::from_stdout(&header, false, bcf::Format::BCF)?
        };
        Ok(self.bcf_writer(writer))
    }
}

impl<L, Pr, Po, ModelPayload> Caller<L, Pr, Po, ModelPayload>
where
    L: bayesian::model::Likelihood<
        ModelPayload,
        Event = AlleleFreqCombination,
        Data = model::modes::generic::Data,
    >,
    Pr: bayesian::model::Prior<Event = AlleleFreqCombination> + model::modes::UniverseDrivenPrior,
    Po: bayesian::model::Posterior<
        BaseEvent = AlleleFreqCombination,
        Event = model::Event,
        Data = model::modes::generic::Data,
    >,
    ModelPayload: Default,
{
    pub fn n_samples(&self) -> usize {
        self.samplenames.len()
    }

    pub fn call(&mut self) -> Result<(), Box<dyn Error>> {
        for obs_reader in self.observations.iter() {
            let mut valid = false;
            for record in obs_reader.header().header_records() {
                if let bcf::HeaderRecord::Generic { key, value } = record {
                    if key == "varlociraptor_observation_format_version" {
                        if value == OBSERVATION_FORMAT_VERSION {
                            valid = true;
                        }
                    }
                }
            }
            if !valid {
                return Err(errors::Error::InvalidObservationFormat)?;
            }
        }

        let mut rid = None;
        let mut events = Vec::new();
        let mut i = 0;
        loop {
            let mut records = self.observations.map(|reader| reader.empty_record());
            let mut eof = Vec::new();
            for (reader, record) in self.observations.iter_mut().zip(records.iter_mut()) {
                eof.push(!reader.read(record)?);
            }

            if eof.iter().all(|v| *v) {
                return Ok(());
            } else if !eof.iter().all(|v| !v) {
                // only some are EOF, this is an error
                return Err(errors::Error::InconsistentObservations)?;
            }

            i += 1;

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
                    return Err(errors::Error::InconsistentObservations)?;
                }
            }

            // rid must not be missind
            let current_rid =
                current_rid.ok_or_else(|| errors::Error::RecordMissingChrom { i: i + 1 })?;

            if !rid.map_or(false, |rid: u32| current_rid == rid) {
                // rid is not the same as before, obtain event universe
                let contig = str::from_utf8(
                    self.observations
                        .first()
                        .unwrap()
                        .header()
                        .rid2name(current_rid)?,
                )
                .unwrap()
                .to_owned();

                // clear old events
                events.clear();

                // register absent event
                events.push(model::Event {
                    name: "absent".to_owned(),
                    vafs: grammar::VAFTree::absent(self.n_samples()),
                    strand_bias: StrandBias::None,
                });

                // add events from scenario
                for (event_name, vaftree) in self.scenario.vaftrees(&contig)? {
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

                self.model.prior_mut().set_universe(vaf_universes.build());

                rid = Some(current_rid);
            }

            let call = self.call_record(&mut records, &events)?;

            if let Some(call) = call {
                call.write_final_record(&mut self.bcf_writer)?;
            }
            if i % 100 == 0 {
                info!("{} records processed.", i);
            }
        }
    }

    fn call_record(
        &mut self,
        records: &mut grammar::SampleInfo<bcf::Record>,
        event_universe: &[Po::Event],
    ) -> Result<Option<Call>, Box<dyn Error>> {
        let (mut call, mut variant_builder, snv) = {
            let first_record = records
                .first_mut()
                .expect("bug: there must be at least one record");
            let start = first_record.pos();
            let chrom = chrom(
                &self
                    .observations
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
                .build()?;

            let mut variant_builder = VariantBuilder::default();
            variant_builder.record(first_record)?;

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

            (call, variant_builder, snv)
        };

        // obtain pileups
        let mut pileups = Vec::new();
        for record in records.iter_mut() {
            pileups.push(read_observations(record)?);
        }
        let data = model::modes::generic::Data::new(pileups, snv);

        // Compute probabilities for given events.
        let m = self.model.compute(event_universe.iter().cloned(), &data);

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
        variant_builder.event_probs(Some(event_probs));

        // add sample specific information
        variant_builder.sample_info(if let Some(map_estimates) = m.maximum_posterior() {
            Some(
                data.into_pileups()
                    .into_iter()
                    .zip(map_estimates.into_iter())
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
        });

        call.variants.push(variant_builder.build()?);

        Ok(Some(call))
    }
}

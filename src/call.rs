// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::HashMap;
use std::error::Error;
use std::path::Path;
use std::str;

use bio::io::fasta;
use bio::stats::{bayesian, LogProb, PHREDProb};
use derive_builder::Builder;
use itertools::join;
use itertools::Itertools;
use rust_htslib::bcf::{self, record::Numeric, Read};
use vec_map::VecMap;

use crate::grammar;
use crate::model;
use crate::model::evidence::observation::expected_depth;
use crate::model::evidence::Observation;
use crate::model::sample::{Pileup, Sample};
use crate::model::{AlleleFreq, StrandBias};
use crate::utils;

pub type AlleleFreqCombination = Vec<model::likelihood::Event>;

#[derive(Default, Clone, Debug, Builder)]
pub struct Call {
    chrom: Vec<u8>,
    pos: u32,
    id: Option<Vec<u8>>,
    #[builder(default = "Vec::new()")]
    variants: Vec<Variant>,
}

#[derive(Default, Clone, Debug, Builder)]
pub struct Variant {
    ref_allele: Vec<u8>,
    alt_allele: Vec<u8>,
    #[builder(default = "None")]
    svlen: Option<i32>,
    event_probs: HashMap<String, LogProb>,
    #[builder(default = "Vec::new()")]
    sample_info: Vec<SampleInfo>,
}

#[derive(Default, Clone, Debug, Builder)]
pub struct SampleInfo {
    allelefreq_estimate: AlleleFreq,
    #[builder(default = "Vec::new()")]
    observations: Vec<Observation>,
    strand_bias: StrandBias,
}

pub fn event_tag_name(event: &str) -> String {
    format!("PROB_{}", event.to_ascii_uppercase())
}

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller<L, Pr, Po, ModelPayload>
where
    L: bayesian::model::Likelihood<ModelPayload>,
    Pr: bayesian::model::Prior,
    Po: bayesian::model::Posterior<Event = model::Event>,
    ModelPayload: Default,
{
    samples: grammar::SampleInfo<Sample>,
    #[builder(private)]
    reference_buffer: utils::ReferenceBuffer,
    #[builder(private)]
    bcf_reader: bcf::Reader,
    #[builder(private)]
    bcf_writer: bcf::Writer,
    #[builder(private)]
    events: HashMap<String, model::Event>,
    #[builder(private)]
    strand_bias_events: Vec<model::Event>,
    model: bayesian::Model<L, Pr, Po, ModelPayload>,
    omit_snvs: bool,
    omit_indels: bool,
    max_indel_len: u32,
}

impl<L, Pr, Po, ModelPayload> CallerBuilder<L, Pr, Po, ModelPayload>
where
    L: bayesian::model::Likelihood<ModelPayload>,
    Pr: bayesian::model::Prior,
    Po: bayesian::model::Posterior<Event = model::Event>,
    ModelPayload: Default,
{
    pub fn reference<P: AsRef<Path>>(self, path: P) -> Result<Self, Box<Error>> {
        Ok(self.reference_buffer(utils::ReferenceBuffer::new(
            fasta::IndexedReader::from_file(&path)?,
        )))
    }

    pub fn inbcf<P: AsRef<Path>>(self, path: Option<P>) -> Result<Self, Box<Error>> {
        Ok(self.bcf_reader(if let Some(path) = path {
            bcf::Reader::from_path(path)?
        } else {
            bcf::Reader::from_stdin()?
        }))
    }

    /// Register events.
    pub fn event(mut self, name: &str, event: grammar::VAFTree) -> Self {
        assert_ne!(
            name.to_ascii_lowercase(),
            "artifact",
            "the event name artifact is reserved for internal use"
        );
        assert_ne!(
            name.to_ascii_lowercase(),
            "absent",
            "the 'absent' event will be created automatically"
        );
        // TODO check that this is not the absent event by looking at !

        if self.events.is_none() {
            let mut events = HashMap::default();
            if let Some(ref samples) = self.samples {
                events.insert(
                    "absent".to_owned(),
                    model::Event {
                        vafs: grammar::VAFTree::absent(samples.len()),
                        strand_bias: StrandBias::None,
                    },
                );
            } else {
                panic!("bug: events must be registered after adding samples");
            }
            self = self.events(events);
            self = self.strand_bias_events(Vec::new());
        }

        let annotate_event = |strand_bias| model::Event {
            vafs: event.clone(),
            strand_bias: strand_bias,
        };

        self.events
            .as_mut()
            .unwrap()
            .insert(name.to_owned(), annotate_event(StrandBias::None));

        // Add strand bias cases.
        self.strand_bias_events
            .as_mut()
            .unwrap()
            .push(annotate_event(StrandBias::Forward));
        self.strand_bias_events
            .as_mut()
            .unwrap()
            .push(annotate_event(StrandBias::Reverse));

        self
    }

    pub fn outbcf<P: AsRef<Path>>(self, path: Option<P>) -> Result<Self, Box<Error>> {
        let mut header = bcf::Header::new();

        // register samples
        for sample in self
            .samples
            .as_ref()
            .expect(".samples() has to be called before .outbcf()")
            .iter()
        {
            header.push_sample(sample.name().as_bytes());
        }

        // register events
        for event in self
            .events
            .as_ref()
            .expect(".events() has to be called before .outbcf()")
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
             Description=\"Posterior probability for strand bias artifact\">",
        );

        // register SVLEN
        header.push_record(
            b"##INFO=<ID=SVLEN,Number=A,Type=Integer,\
              Description=\"Difference in length between REF and ALT alleles\">",
        );

        // register sample specific tags
        header.push_record(
            b"##FORMAT=<ID=DP,Number=A,Type=Integer,\
              Description=\"Expected sequencing depth, while considering mapping uncertainty\">",
        );
        header.push_record(
            b"##FORMAT=<ID=AF,Number=A,Type=Float,\
              Description=\"Maximum a posteriori probability estimate of allele frequency\">",
        );
        header.push_record(
            b"##FORMAT=<ID=OBS,Number=A,Type=String,\
              Description=\"Posterior odds for alt allele of each fragment as Kass Raftery \
              scores: N=none, B=barely, P=positive, S=strong, V=very strong (lower case if \
              probability for correct mapping of fragment is <95%)\">",
        );
        header.push_record(
            b"##FORMAT=<ID=SB,Number=A,Type=String,\
              Description=\"Strand bias estimate: + indicates that ALT allele is associated with \
              forward strand, - indicates that ALT allele is associated with reverse strand, \
              - indicates no strand bias.\">",
        );

        // register sequences
        for sequence in self
            .reference_buffer
            .as_ref()
            .expect(".reference() has to be called before .outbcf()")
            .reader
            .index
            .sequences()
        {
            header.push_record(
                format!("##contig=<ID={},length={}>", sequence.name, sequence.len).as_bytes(),
            );
        }

        let writer = if let Some(path) = path {
            bcf::Writer::from_path(path, &header, false, false)?
        } else {
            bcf::Writer::from_stdout(&header, false, false)?
        };
        Ok(self.bcf_writer(writer))
    }
}

impl<L, Pr, Po, ModelPayload> Caller<L, Pr, Po, ModelPayload>
where
    L: bayesian::model::Likelihood<ModelPayload, Event = AlleleFreqCombination, Data = Vec<Pileup>>,
    Pr: bayesian::model::Prior<Event = AlleleFreqCombination>,
    Po: bayesian::model::Posterior<
        BaseEvent = AlleleFreqCombination,
        Event = model::Event,
        Data = Vec<Pileup>,
    >,
    ModelPayload: Default,
{
    pub fn call(&mut self) -> Result<(), Box<Error>> {
        let mut universe = self.events.values().cloned().collect_vec();
        universe.extend(self.strand_bias_events.iter().cloned());

        let mut i = 0;
        let mut record = self.bcf_reader.empty_record();
        loop {
            match self.bcf_reader.read(&mut record) {
                Err(bcf::ReadError::NoMoreRecord) => return Ok(()),
                Err(e) => return Err(Box::new(e)),
                Ok(()) => (),
            }

            i += 1;
            let call = self.call_record(&mut record, &universe)?;
            if let Some(call) = call {
                self.write_call(&call)?;
            }
            if i % 100 == 0 {
                info!("{} records processed.", i);
            }
        }
    }

    fn write_call(&mut self, call: &Call) -> Result<(), Box<Error>> {
        let rid = self.bcf_writer.header().name2rid(&call.chrom)?;
        for (ref_allele, group) in call
            .variants
            .iter()
            .group_by(|variant| &variant.ref_allele)
            .into_iter()
        {
            let mut record = self.bcf_writer.empty_record();
            record.set_rid(&Some(rid));
            record.set_pos(call.pos as i32);
            // set ID if present
            if let Some(ref id) = call.id {
                record.set_id(id)?;
            }

            let mut event_probs = HashMap::new();
            let mut allelefreq_estimates = VecMap::new();
            let mut observations = VecMap::new();
            let mut obs_counts = VecMap::new();
            let mut strand_bias = VecMap::new();
            let mut alleles = Vec::new();
            let mut svlens = Vec::new();
            alleles.push(&ref_allele[..]);

            // collect per group information
            for variant in group {
                alleles.push(&variant.alt_allele[..]);

                for (event, prob) in &variant.event_probs {
                    event_probs
                        .entry(event)
                        .or_insert_with(|| Vec::new())
                        .push(*prob);
                }

                for (i, sample_info) in variant.sample_info.iter().enumerate() {
                    strand_bias.entry(i).or_insert_with(Vec::new).push(
                        match sample_info.strand_bias {
                            StrandBias::None => '.',
                            StrandBias::Forward => '+',
                            StrandBias::Reverse => '-',
                        },
                    );

                    allelefreq_estimates
                        .entry(i)
                        .or_insert_with(|| Vec::new())
                        .push(*sample_info.allelefreq_estimate as f32);

                    obs_counts
                        .entry(i)
                        .or_insert_with(|| Vec::new())
                        .push(expected_depth(&sample_info.observations) as i32);

                    observations.entry(i).or_insert_with(|| Vec::new()).push({
                        utils::generalized_cigar(
                            sample_info.observations.iter().map(|obs| {
                                let score = utils::evidence_kass_raftery_to_letter(
                                    obs.bayes_factor_alt().evidence_kass_raftery(),
                                );
                                format!(
                                    "{}{}",
                                    if obs.prob_mapping < LogProb(0.95_f64.ln()) {
                                        score.to_ascii_lowercase()
                                    } else {
                                        score.to_ascii_uppercase()
                                    },
                                    match (obs.forward_strand, obs.reverse_strand) {
                                        (true, true) => '*',
                                        (false, true) => '-',
                                        (true, false) => '+',
                                        _ => panic!("bug: unknown strandedness"),
                                    }
                                )
                            }),
                            false,
                        )
                    })
                }

                if let Some(svlen) = variant.svlen {
                    svlens.push(svlen);
                } else {
                    svlens.push(i32::missing());
                }
            }

            // set alleles
            record.set_alleles(&alleles)?;
            record.push_info_integer(b"SVLEN", &svlens)?;

            // set qual
            record.set_qual(f32::missing());

            // set event probabilities
            for (event, probs) in event_probs {
                let probs = probs
                    .iter()
                    .map(|p| PHREDProb::from(*p).abs() as f32)
                    .collect_vec();
                record.push_info_float(event_tag_name(event).as_bytes(), &probs)?;
            }

            // set sample info
            let dp = obs_counts.values().flatten().cloned().collect_vec();
            record.push_format_integer(b"DP", &dp)?;

            let afs = allelefreq_estimates
                .values()
                .flatten()
                .cloned()
                .collect_vec();
            record.push_format_float(b"AF", &afs)?;

            let obs = if observations.values().any(|obs| !obs.is_empty()) {
                observations
                    .values()
                    .map(|obs| join(obs.iter(), ",").into_bytes())
                    .collect_vec()
            } else {
                vec![b".".to_vec()]
            };
            record.push_format_string(b"OBS", &obs)?;

            let sb = strand_bias
                .values()
                .map(|sb| join(sb.iter(), ",").into_bytes())
                .collect_vec();

            record.push_format_string(b"SB", &sb)?;

            self.bcf_writer.write(&record)?;
        }

        Ok(())
    }

    fn call_record(
        &mut self,
        record: &mut bcf::Record,
        universe: &[Po::Event],
    ) -> Result<Option<Call>, Box<Error>> {
        let start = record.pos();
        let chrom = chrom(&self.bcf_reader, &record);
        let variants = utils::collect_variants(
            record,
            self.omit_snvs,
            self.omit_indels,
            Some(0..self.max_indel_len + 1),
        )?;

        if variants.is_empty() || variants.iter().all(|v| v.is_none()) {
            return Ok(None);
        }

        let mut call = CallBuilder::default()
            .chrom(chrom.to_owned())
            .pos(start)
            .id({
                let id = record.id();
                if id == b"." {
                    None
                } else {
                    Some(id)
                }
            })
            .variants(Vec::new())
            .build()?;

        for variant in variants.into_iter() {
            if let Some(variant) = variant {
                let mut pileups = Vec::new();
                for sample in self.samples.iter_mut() {
                    let chrom_seq = self.reference_buffer.seq(&chrom)?;
                    let pileup = sample.extract_observations(start, &variant, chrom, chrom_seq)?;
                    pileups.push(pileup);
                }

                // Compute probabilities for given events.
                let m = self.model.compute(universe.iter().cloned(), &pileups);
                let mut variant_builder = VariantBuilder::default();

                // add calling results
                let mut event_probs: HashMap<String, LogProb> = self
                    .events
                    .iter()
                    .map(|(name, event)| (name.clone(), m.posterior(event).unwrap()))
                    .collect();
                // generate artifact event
                event_probs.insert(
                    "artifact".to_owned(),
                    LogProb::ln_sum_exp(
                        &self
                            .strand_bias_events
                            .iter()
                            .map(|event| m.posterior(event).unwrap())
                            .collect_vec(),
                    ),
                );
                variant_builder.event_probs(event_probs);

                // add sample specific information
                variant_builder.sample_info(if let Some(map_estimates) = m.maximum_posterior() {
                    pileups
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
                        .collect_vec()
                } else {
                    // no observations
                    Vec::new()
                });

                let start = start as usize;

                // add variant information
                let chrom_seq = self.reference_buffer.seq(&chrom)?;
                match variant {
                    model::Variant::Deletion(l) => {
                        let svlen = -(l as i32);
                        if l <= 50 {
                            variant_builder
                                .ref_allele(
                                    chrom_seq[start..start + 1 + l as usize].to_ascii_uppercase(),
                                )
                                .alt_allele(chrom_seq[start..start + 1].to_ascii_uppercase())
                                .svlen(Some(svlen));
                        } else {
                            variant_builder
                                .ref_allele(chrom_seq[start..start + 1].to_ascii_uppercase())
                                .alt_allele(b"<DEL>".to_ascii_uppercase())
                                .svlen(Some(svlen));
                        }
                    }
                    model::Variant::Insertion(ref seq) => {
                        let svlen = seq.len() as i32;
                        let ref_allele = vec![chrom_seq[start]];
                        let mut alt_allele = ref_allele.clone();
                        alt_allele.extend(seq);

                        variant_builder
                            .ref_allele(ref_allele.to_ascii_uppercase())
                            .alt_allele(alt_allele.to_ascii_uppercase())
                            .svlen(Some(svlen));
                    }
                    model::Variant::SNV(base) => {
                        variant_builder
                            .ref_allele(chrom_seq[start..start + 1].to_ascii_uppercase())
                            .alt_allele(vec![base].to_ascii_uppercase());
                    }
                    model::Variant::None => {
                        variant_builder
                            .ref_allele(chrom_seq[start..start + 1].to_ascii_uppercase())
                            .alt_allele(b"<REF>".to_ascii_uppercase());
                    }
                }

                call.variants.push(variant_builder.build()?);
            }
        }

        Ok(Some(call))
    }
}

fn chrom<'a>(inbcf: &'a bcf::Reader, record: &bcf::Record) -> &'a [u8] {
    inbcf.header().rid2name(record.rid().unwrap())
}

// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::HashMap;
use std::error::Error;
use std::path::Path;
use std::str;

use bio::io::fasta;
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::{bayesian, LogProb, PHREDProb};
use counter::Counter;
use derive_builder::Builder;
use itertools::Itertools;
use rust_htslib::bcf::{self, record::Numeric, Read};
use vec_map::VecMap;

use crate::model;
use crate::model::evidence::Observation;
use crate::model::sample::{Pileup, Sample};
use crate::model::{AlleleFreq, AlleleFreqs};
use crate::utils;

pub type AlleleFreqCombination = Vec<AlleleFreq>;

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
}

impl SampleInfo {
    fn observation_count(&self, alt_allele: bool) -> u32 {
        self.observations
            .iter()
            .map(|obs| {
                let bf = if alt_allele {
                    obs.bayes_factor_alt()
                } else {
                    obs.bayes_factor_ref()
                };
                if bf.evidence_kass_raftery() >= KassRaftery::Positive {
                    1
                } else {
                    0
                }
            })
            .sum()
    }

    pub fn n_obs_ref(&self) -> u32 {
        self.observation_count(false)
    }

    pub fn n_obs_alt(&self) -> u32 {
        self.observation_count(true)
    }
}

pub fn event_tag_name(event: &str) -> String {
    format!("PROB_{}", event.to_ascii_uppercase())
}

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller<A, L, Pr, Po>
where
    A: AlleleFreqs,
    L: bayesian::model::Likelihood,
    Pr: bayesian::model::Prior,
    Po: bayesian::model::Posterior,
{
    samples: Vec<Sample>,
    #[builder(private)]
    reference_buffer: utils::ReferenceBuffer,
    #[builder(private)]
    bcf_reader: bcf::Reader,
    #[builder(private)]
    bcf_writer: bcf::Writer,
    #[builder(private)]
    events: HashMap<String, Vec<A>>,
    model: bayesian::Model<L, Pr, Po>,
    omit_snvs: bool,
    omit_indels: bool,
    max_indel_len: u32,
    exclusive_end: bool,
}

impl<A, L, Pr, Po> CallerBuilder<A, L, Pr, Po>
where
    A: AlleleFreqs,
    L: bayesian::model::Likelihood,
    Pr: bayesian::model::Prior,
    Po: bayesian::model::Posterior<Event = Vec<A>>,
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

    pub fn event<E: Into<Vec<A>>>(mut self, name: &str, event: E) -> Self {
        if self.events.is_none() {
            self = self.events(HashMap::default());
        }
        self.events
            .as_mut()
            .unwrap()
            .insert(name.to_owned(), event.into());

        self
    }

    pub fn outbcf<P: AsRef<Path>>(self, path: Option<P>) -> Result<Self, Box<Error>> {
        let mut header = bcf::Header::new();

        // register samples
        for sample in self
            .samples
            .as_ref()
            .expect(".samples() has to be called before .outbcf()")
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
                     Description=\"Posterior probability for event {}\">",
                    event_tag_name(event),
                    event
                )
                .as_bytes(),
            );
        }

        // register SVLEN
        header.push_record(
            b"##INFO=<ID=SVLEN,Number=A,Type=Integer,\
              Description=\"Difference in length between REF and ALT alleles\">",
        );

        // register sample specific tags
        header.push_record(
            b"##FORMAT=<ID=AF,Number=A,Type=Float,\
              Description=\"Maximum a posteriori probability estimate of allele frequency\">",
        );
        header.push_record(
            b"##FORMAT=<ID=OBS,Number=A,Type=String,\
              Description=\"Posterior odds for alt allele of each fragment as Kass Raftery scores: N=none, B=barely, P=positive, S=strong, V=very strong (lower case if probability for correct mapping of fragment is <95%)\">",
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

impl<A, L, Pr, Po> Caller<A, L, Pr, Po>
where
    A: AlleleFreqs + Clone + Ord,
    L: bayesian::model::Likelihood<Event = AlleleFreqCombination, Data = Vec<Pileup>>,
    Pr: bayesian::model::Prior<Event = AlleleFreqCombination>,
    Po: bayesian::model::Posterior<
        BaseEvent = AlleleFreqCombination,
        Event = Vec<A>,
        Data = Vec<Pileup>,
    >,
{
    pub fn call(&mut self) -> Result<(), Box<Error>> {
        let mut i = 0;
        loop {
            let record = self.next_record()?;
            if let Some(mut record) = record {
                let call = self.call_record(&mut record)?;
                if let Some(call) = call {
                    self.write_call(&call)?;
                    i += 1;
                    if i % 100 == 0 {
                        info!("{} records written.", i);
                    }
                }
            } else {
                // done
                return Ok(());
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
                    allelefreq_estimates
                        .entry(i)
                        .or_insert_with(|| Vec::new())
                        .push(*sample_info.allelefreq_estimate as f32);
                    observations.entry(i).or_insert_with(|| Vec::new()).push({
                        let obs: Counter<u8> = sample_info
                            .observations
                            .iter()
                            .map(|obs| {
                                let score = match obs.bayes_factor_alt().evidence_kass_raftery() {
                                    KassRaftery::Barely => b'B',
                                    KassRaftery::None => b'N',
                                    KassRaftery::Positive => b'P',
                                    KassRaftery::Strong => b'S',
                                    KassRaftery::VeryStrong => b'V',
                                };
                                if obs.prob_mapping < LogProb(0.95_f64.ln()) {
                                    score.to_ascii_lowercase()
                                } else {
                                    score.to_ascii_uppercase()
                                }
                            })
                            .collect();

                        obs.most_common().into_iter().map(|(score, count)| format!("{}{}", count, score as char).into_bytes()).flatten().collect_vec()
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
            let afs = allelefreq_estimates
                .values()
                .cloned()
                .flatten()
                .collect_vec();
            record.push_format_float(b"AF", &afs)?;

            let obs = observations
                .values()
                .flatten()
                .map(|obs| &obs[..])
                .collect_vec();
            record.push_format_string(b"OBS", &obs)?;

            self.bcf_writer.write(&record)?;
        }

        Ok(())
    }

    fn next_record(&mut self) -> Result<Option<bcf::Record>, Box<Error>> {
        let mut record = self.bcf_reader.empty_record();
        match self.bcf_reader.read(&mut record) {
            Err(bcf::ReadError::NoMoreRecord) => Ok(None),
            Err(e) => Err(Box::new(e)),
            Ok(()) => Ok(Some(record)),
        }
    }

    fn call_record(&mut self, record: &mut bcf::Record) -> Result<Option<Call>, Box<Error>> {
        let start = record.pos();
        let chrom = chrom(&self.bcf_reader, &record);
        let variants = utils::collect_variants(
            record,
            self.omit_snvs,
            self.omit_indels,
            Some(0..self.max_indel_len),
            self.exclusive_end,
        )?;

        if variants.is_empty() {
            return Ok(None);
        }

        let chrom_seq = self.reference_buffer.seq(&chrom)?;

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
                for sample in &mut self.samples {
                    let pileup = sample.extract_observations(start, &variant, chrom, chrom_seq)?;
                    pileups.push(pileup);
                }

                // compute probabilities
                let m = self.model.compute(self.events.values().cloned(), &pileups);

                let mut variant_builder = VariantBuilder::default();

                // add calling results
                variant_builder.event_probs(
                    self.events
                        .iter()
                        .map(|(name, event)| (name.clone(), m.posterior(event).unwrap()))
                        .collect(),
                );

                // add sample specific information
                variant_builder.sample_info(if let Some(map_estimates) = m.maximum_posterior() {
                    pileups
                        .into_iter()
                        .zip(map_estimates.clone().into_iter())
                        .map(|(pileup, estimate)| {
                            SampleInfoBuilder::default()
                                .allelefreq_estimate(estimate.clone())
                                .observations(pileup)
                                .build()
                                .unwrap()
                        })
                        .collect_vec()
                } else {
                    // no observations
                    Vec::new()
                });

                let start = start as usize;

                // add variant information
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

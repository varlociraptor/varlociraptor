use std::cmp::Ord;
use std::error::Error;
use std::path::Path;

use bio::io::fasta;
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::{bayesian, LogProb, PHREDProb};
use csv;
use derive_builder::Builder;
use itertools::Itertools;
use rust_htslib::bcf::{self, Read};
use vec_map::VecMap;

use crate::model;
use crate::model::evidence::Observation;
use crate::model::sample::{Pileup, Sample};
use crate::model::AlleleFreq;
use crate::utils;
use crate::Event;

pub struct PairEvent<A: AlleleFreqs, B: AlleleFreqs> {
    /// event name
    pub name: String,
    /// allele frequencies for case sample
    pub af_case: A,
    /// allele frequencies for control sample
    pub af_control: B,
}

impl<A: AlleleFreqs, B: AlleleFreqs> Event for PairEvent<A, B> {
    fn name(&self) -> &str {
        &self.name
    }
}

#[derive(Default, Clone, Debug, Builder)]
pub struct Call {
    chrom: Vec<u8>,
    pos: u32,
    variants: Vec<Variant>,
}

#[derive(Default, Clone, Debug, Builder)]
pub struct Variant {
    ref_allele: Vec<u8>,
    alt_allele: Vec<u8>,
    svlen: Option<u32>,
    event_probs: Vec<LogProb>,
    sample_info: Vec<SampleInfo>,
}

#[derive(Default, Clone, Debug, Builder)]
pub struct SampleInfo {
    allelefreq_estimate: AlleleFreq,
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

pub struct BCFWriter<E: Event + Clone> {
    bcf: bcf::Writer,
    events: Vec<E>,
}

impl<E: Event + Clone> BCFWriter<E> {
    pub fn new<P: AsRef<Path>>(
        path: Option<P>,
        sample_names: &[Vec<u8>],
        events: &[E],
        sequences: &[fasta::Sequence],
    ) -> Result<Self, Box<Error>> {
        let mut header = bcf::Header::new();
        for sample in sample_names {
            header.push_sample(sample);
        }
        for event in events {
            header.push_record(
                event
                    .header_entry("PROB", "Posterior probability for ")
                    .as_bytes(),
            );
        }
        header.push_record(
            b"##FORMAT=<ID=AF,Number=A,Type=Float,\
            Description=\"Maximum a posteriori probability estimate of allele frequency.\">",
        );
        for sequence in sequences {
            header.push_record(
                format!("##contig=<ID={},length={}>", sequence.name, sequence.len).as_bytes(),
            );
        }

        let bcf = if let Some(path) = path {
            bcf::Writer::from_path(path, &header, false, false)?
        } else {
            bcf::Writer::from_stdout(&header, false, false)?
        };
        let events = events.to_vec();

        Ok(BCFWriter { bcf, events })
    }

    pub fn write(&mut self, call: &Call) -> Result<(), Box<Error>> {
        let rid = self.bcf.header().name2rid(&call.chrom)?;
        for (ref_allele, group) in call
            .variants
            .iter()
            .group_by(|variant| &variant.ref_allele)
            .into_iter()
        {
            let mut record = self.bcf.empty_record();
            record.set_rid(&Some(rid));
            record.set_pos(call.pos as i32);

            let mut event_probs = vec![Vec::new(); self.events.len()];
            let mut allelefreq_estimates = VecMap::new();
            let mut alleles = Vec::new();
            alleles.push(&ref_allele[..]);

            // collect per group information
            for variant in group {
                alleles.push(&variant.alt_allele[..]);

                for i in 0..self.events.len() {
                    event_probs[i].push(variant.event_probs[i]);
                }
                for (i, sample_info) in variant.sample_info.iter().enumerate() {
                    allelefreq_estimates
                        .entry(i)
                        .or_insert_with(|| Vec::new())
                        .push(*sample_info.allelefreq_estimate as f32);
                }
            }

            // set alleles
            record.set_alleles(&alleles)?;
            // set event probabilities
            for (event, probs) in self.events.iter().zip(event_probs) {
                let probs = probs
                    .iter()
                    .map(|p| PHREDProb::from(*p).abs() as f32)
                    .collect_vec();
                record.push_info_float(&event.tag_name("PROB").as_bytes(), &probs)?;
            }
            // set sample info
            let afs = allelefreq_estimates
                .values()
                .cloned()
                .flatten()
                .collect_vec();
            record.push_format_float(b"AF", &afs)?;

            self.bcf.write(&record)?;
        }

        Ok(())
    }
}

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller<L, Pr, Po>
where
    L: bayesian::model::Likelihood,
    Pr: bayesian::model::Prior,
    Po: bayesian::model::Posterior,
{
    samples: Vec<Sample>,
    reference_buffer: utils::ReferenceBuffer,
    candidates: bcf::Reader,
    events: Vec<Po::Event>,
    model: bayesian::Model<L, Pr, Po>,
    omit_snvs: bool,
    omit_indels: bool,
    max_indel_len: Option<u32>,
    exclusive_end: bool,
}

impl<AlleleFreqCombination, Event, L, Pr, Po> Caller<L, Pr, Po>
where
    AlleleFreqCombination: Ord + Clone + IntoIterator<Item = AlleleFreq>,
    Event: Ord + Clone,
    L: bayesian::model::Likelihood<Event = AlleleFreqCombination, Data = Vec<Pileup>>,
    Pr: bayesian::model::Prior<Event = AlleleFreqCombination>,
    Po: bayesian::model::Posterior<
        BaseEvent = AlleleFreqCombination,
        Event = Event,
        Data = Vec<Pileup>,
    >,
{
    fn call(&mut self) -> Result<Option<Call>, Box<Error>> {
        let mut record = self.candidates.empty_record();
        match self.candidates.read(&mut record) {
            Err(bcf::ReadError::NoMoreRecord) => return Ok(None),
            Err(e) => return Err(Box::new(e)),
            Ok(()) => (),
        }

        let start = record.pos();
        let chrom = chrom(&self.candidates, &record);
        let variants = utils::collect_variants(
            &mut record,
            self.omit_snvs,
            self.omit_indels,
            self.max_indel_len.map(|l| 0..l),
            self.exclusive_end,
        )?;

        let chrom_seq = self.reference_buffer.seq(&chrom)?;

        let mut call = CallBuilder::default()
            .chrom(chrom.to_owned())
            .pos(start)
            .build()?;

        for variant in variants.into_iter() {
            if let Some(variant) = variant {
                let mut pileups = Vec::new();
                for sample in &mut self.samples {
                    let pileup = sample.extract_observations(start, &variant, chrom, chrom_seq)?;
                    pileups.push(pileup);
                }

                // compute probabilities
                // TODO remove type annotation
                let m = self.model.compute(&self.events, &pileups);

                let mut variant_builder = VariantBuilder::default();

                // add calling results
                variant_builder.event_probs(
                    self.events
                        .iter()
                        .map(|event| m.posterior(event).unwrap())
                        .collect_vec(),
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
                        if l <= 10 {
                            variant_builder
                                .ref_allele(chrom_seq[start - 1..start + l as usize].to_vec())
                                .alt_allele(chrom_seq[start - 1..start].to_vec());
                        } else {
                            variant_builder
                                .ref_allele(chrom_seq[start..start + 1].to_vec())
                                .alt_allele(b"<DEL>".to_vec())
                                .svlen(Some(l));
                        }
                    }
                    model::Variant::Insertion(ref seq) => {
                        let ref_allele = vec![chrom_seq[start - 1]];
                        let mut alt_allele = ref_allele.clone();
                        alt_allele.extend(seq);

                        variant_builder
                            .ref_allele(ref_allele)
                            .alt_allele(alt_allele);
                    }
                    model::Variant::SNV(base) => {
                        variant_builder
                            .ref_allele(chrom_seq[start..start + 1].to_vec())
                            .alt_allele(vec![base]);
                    }
                    model::Variant::None => {
                        variant_builder
                            .ref_allele(chrom_seq[start..start + 1].to_vec())
                            .alt_allele(b"<REF>".to_vec());
                    }
                }

                call.variants.push(variant_builder.build()?);
            }
        }

        Ok(Some(call))
    }
}

impl<AlleleFreqCombination, Event, L, Pr, Po> Iterator for Caller<L, Pr, Po>
where
    AlleleFreqCombination: Ord + Clone + IntoIterator<Item = AlleleFreq>,
    Event: Ord + Clone,
    L: bayesian::model::Likelihood<Event = AlleleFreqCombination, Data = Vec<Pileup>>,
    Pr: bayesian::model::Prior<Event = AlleleFreqCombination>,
    Po: bayesian::model::Posterior<
        BaseEvent = AlleleFreqCombination,
        Event = Event,
        Data = Vec<Pileup>,
    >,
{
    type Item = Result<Call, Box<Error>>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.call() {
            Ok(Some(call)) => Some(Ok(call)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

fn chrom<'a>(inbcf: &'a bcf::Reader, record: &bcf::Record) -> &'a [u8] {
    inbcf.header().rid2name(record.rid().unwrap())
}

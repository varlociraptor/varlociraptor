use std::cmp::Ord;
use std::error::Error;

use bio::stats::{bayesian, LogProb};
use csv;
use derive_builder::Builder;
use itertools::Itertools;
use rust_htslib::bcf::{self, Read};

use crate::model;
use crate::model::sample::{Pileup, Sample};
use crate::model::evidence::Observation;
use crate::model::AlleleFreq;
use crate::utils;

#[derive(Default, Clone, Debug, Builder)]
pub struct Call {
    chrom: Vec<u8>,
    pos: usize,
    variants: Vec<Variant>,
}

#[derive(Default, Clone, Debug, Builder)]
pub struct Variant {
    ref_allele: Vec<u8>,
    alt_allele: Vec<u8>,
    svlen: Option<usize>,
    event_probs: Vec<LogProb>,
    sample_info: Vec<SampleInfo>,
}

#[derive(Default, Clone, Debug, Builder)]
pub struct SampleInfo {
    maximum_posterior_estimate: AlleleFreq,
    observations: Vec<Observation>,
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
    inbcf: bcf::Reader,
    events: Vec<Po::Event>,
    model: bayesian::Model<L, Pr, Po>,
    omit_snvs: bool,
    omit_indels: bool,
    max_indel_len: Option<u32>,
    exclusive_end: bool,
}

impl<AlleleFreqCombination, Event, L, Pr, Po> Caller<L, Pr, Po>
where
    AlleleFreqCombination: Ord + Clone,
    Event: Ord + Clone,
    L: bayesian::model::Likelihood<Event=AlleleFreqCombination, Data=Vec<Pileup>>,
    Pr: bayesian::model::Prior<Event=AlleleFreqCombination>,
    Po: bayesian::model::Posterior<BaseEvent=AlleleFreqCombination, Event=Event, Data=Vec<Pileup>>,

{
    fn call(&mut self) -> Result<Option<Call>, Box<Error>> {
        let mut record = self.inbcf.empty_record();
        match self.inbcf.read(&mut record) {
            Err(bcf::ReadError::NoMoreRecord) => return Ok(None),
            Err(e) => return Err(Box::new(e)),
            Ok(()) => ()
        }

        let start = record.pos();
        let chrom = chrom(&self.inbcf, &record);
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
            .pos(start as usize)
            .build()?;

        for (i, variant) in variants.into_iter().enumerate() {
            if let Some(variant) = variant {
                let pileups = Vec::new();
                for sample in &mut self.samples {
                    let pileup = sample.extract_observations(start, &variant, chrom, chrom_seq)?;
                    pileups.push(pileup);
                }

                // compute probabilities
                let m = self.model.compute(&self.events, &pileups);

                let mut variant_builder = VariantBuilder::default();

                // add calling results
                variant_builder.event_probs(
                    self.events.iter().map(|event| m.posterior(event).unwrap()).collect_vec()
                );

                // add sample specific information
                variant_builder.sample_info(
                    pileups.iter().zip(m.maximum_posterior()).map(|(pileup, estimate)| {
                        SampleInfo::default()
                            .maximum_posterior_estimate(estimate)
                            .observations(pileup)
                    }).collect_vec()
                );

                let start = start as usize;

                // add variant information
                match variant {
                    model::Variant::Deletion(l) => {
                        variant_builder
                            .ref_allele(chrom_seq[start].to_owned())
                            .alt_allele(b"<DEL>".to_owned())
                            .svlen(Some(l));
                    }
                    model::Variant::Insertion(ref seq) => {
                        variant_builder
                            .ref_allele(chrom_seq[start - 1].to_owned())
                            .alt_allele(seq.to_owned())
                            .svlen(Some(seq.len()));
                    }
                    model::Variant::SNV(base) => {
                        variant_builder
                            .ref_allele(chrom_seq[start].to_owned())
                            .alt_allele(vec![base]);
                    }
                }

                call.variants.push(variant_builder.build()?);
            }
        }

        Ok(Some(call))
    }
}


impl<L, Pr, Po> Iterator for Caller<L, Pr, Po>
where
    L: bayesian::model::Likelihood,
    Pr: bayesian::model::Prior,
    Po: bayesian::model::Posterior,
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

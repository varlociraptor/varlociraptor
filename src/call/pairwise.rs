// TODO deprecated, remove.

use std::error::Error;
use std::f32;
use std::path::Path;
use std::str;
use std::cmp::Ordering;

use bio::io::fasta;
use bio::stats::{LogProb, PHREDProb, bayesian};
use csv;
use itertools::Itertools;
use ndarray::prelude::*;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read};
use derive_builder::Builder;

use crate::model::{self, sample::{Sample, Pileup}};
use crate::model::{AlleleFreqs, AlleleFreq};
use crate::model::PairCaller;
use crate::utils;
use crate::Event;
use super::{Call, CallBuilder, VariantBuilder};


#[derive(Builder)]
pub struct Caller<L, Pr, Po>
where
    Po: bayesian::model::Posterior
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
    bcf_record: bcf::Record,
}


impl<L, Pr, Po> Iterator for Caller<L, Pr, Po> {
    type Item = Result<Call, Box<Error>>;

    fn next(&mut self) -> Self::Item {
        self.inbcf.read(self.bcf_record);

        let start = record.pos();
        let chrom = chrom(&self.inbcf, self.bcf_record);
        let variants = utils::collect_variants(
            self.bcf_record,
            self.omit_snvs,
            self.omit_indels,
            self.max_indel_len.map(|l| 0..l),
            self.exclusive_end
        )?;

        let chrom_seq = self.reference_buffer.seq(&chrom)?;

        let pileup = |start, variant, sample| {
            sample.extract_observations(
                start,
                variant,
                chrom,
                chrom_seq,
            )
        };

        let mut call = CallBuilder::default()
            .chrom(chrom)
            .pos(start as usize)
            .build()?;

        for (i, variant) in variants.enumerate() {
            if let Some(variant) = variant {

                let pileups = self.samples.iter_mut().map(|sample|
                    sample.extract_observations(
                        start,
                        variant,
                        chrom,
                        chrom_seq,
                    )
                );
                // TODO remove this constraint
                assert_eq!(pileups.len(), 2, "calling only supported for sample pairs right now.");
                let pileups = (pileups[0], pileups[1]);

                // compute probabilities
                let m = self.model.compute(&self.events, &pileups);

                let mut variant_builder = VariantBuilder::default();

                // add calling results
                variant_build
                    .event_probs(self.events.iter().map(|event| m.posterior(event).unwrap()));
                    .maximum_posterior_estimates(m.maximum_posterior_estimates().to_vec());

                // add variant information
                match variant {
                    Variant::Deletion(l) => {
                        variant_builder
                            .ref_allele(chrom_seq[start])
                            .alt_allele(b"<DEL>")
                            .svlen(Some(l));
                    },
                    Variant::Insertion(ref seq) => {
                        variant_builder
                            .ref_allele(chrom_seq[start - 1])
                            .alt_allele(seq)
                            .svlen(Some(seq.len()));
                    },
                    Variant::SNV(base) => {
                        variant_builder
                            .ref_allele(chrom_seq[start])
                            .alt_allele(vec![base]);
                    }
                }
                variant = variant_build.build()?;

                call.variants.push(variant);
            }
        }

        Ok(call)
    }
}



fn phred_scale<'a, I: IntoIterator<Item = &'a Option<LogProb>>>(probs: I) -> Vec<f32> {
    probs
        .into_iter()
        .map(|&p| match p {
            Some(p) => PHREDProb::from(p).abs() as f32,
            None => f32::missing(),
        })
        .collect_vec()
}

fn pileups<'a, A, B, P>(
    inbcf: &bcf::Reader,
    record: &mut bcf::Record,
    joint_model: &'a mut PairCaller<A, B, P>,
    reference_buffer: &mut utils::ReferenceBuffer,
    omit_snvs: bool,
    omit_indels: bool,
    max_indel_len: Option<u32>,
    exclusive_end: bool,
) -> Result<Vec<Option<model::PairPileup<'a, A, B, P>>>, Box<Error>>
where
    A: AlleleFreqs,
    B: AlleleFreqs,
    P: priors::PairModel<A, B>,
{
    let chrom = chrom(&inbcf, &record);
    let variants = utils::collect_variants(
        record,
        omit_snvs,
        omit_indels,
        max_indel_len.map(|l| 0..l),
        exclusive_end
    )?;

    let chrom_seq = reference_buffer.seq(&chrom)?;

    let mut pileups = Vec::with_capacity(variants.len());
    for variant in variants {
        pileups.push(if let Some(variant) = variant {
            Some(joint_model.pileup(
                chrom,
                record.pos(),
                variant,
                chrom_seq
            )?)
        } else {
            None
        });
    }

    Ok(pileups)
}

/// Call variants with the given model.
///
/// # Arguments
///
/// * `inbcf` - path to BCF/VCF with preprocessed variant calls (None for STDIN).
/// * `outbcf` - path to BCF/VCF with results (None for STDOUT).
/// * `events` - events to call (these have to cover the entire event space!!)
/// * `joint_model` - calling model to use
/// * `omit_snvs` - omit single nucleotide variants
/// * `omit_indels` - omit indels
/// * `outobs` - optional path where to store observations as JSON
///
/// # Returns
///
/// `Result` object with eventual error message.

pub fn call<A, B, P, M, R, W, X, F>(
    inbcf: Option<R>,
    outbcf: Option<W>,
    fasta: &F,
    events: &[PairEvent<A, B>],
    pair_model: &mut PairCaller<A, B, P>,
    omit_snvs: bool,
    omit_indels: bool,
    max_indel_len: Option<u32>,
    outobs: Option<&X>,
    exclusive_end: bool,
) -> Result<(), Box<Error>>
where
    A: AlleleFreqs,
    B: AlleleFreqs,
    P: priors::PairModel<A, B>,
    R: AsRef<Path>,
    W: AsRef<Path>,
    X: AsRef<Path>,
    F: AsRef<Path>,
{
    let fasta = fasta::IndexedReader::from_file(fasta)?;
    let mut reference_buffer = utils::ReferenceBuffer::new(fasta);

    let mut inbcf = match inbcf {
        Some(f) => bcf::Reader::from_path(f)?,
        None => bcf::Reader::from_stdin()?,
    };

    let mut header = bcf::Header::from_template(inbcf.header());
    for event in events {
        header.push_record(
            event
                .header_entry("PROB", "PHRED-scaled probability for")
                .as_bytes(),
        );
    }

    // add tags for expected allele frequency
    header.push_record(
        b"##INFO=<ID=CASE_AF,Number=A,Type=Float,\
        Description=\"Maximum a posteriori probability estimate of allele frequency in case sample.\">"
    );
    header.push_record(
        b"##INFO=<ID=CONTROL_AF,Number=A,Type=Float,\
        Description=\"Maximum a posteriori probability estimate of allele frequency in control sample.\">"
    );

    let mut outbcf = match outbcf {
        Some(f) => bcf::Writer::from_path(f, &header, false, false)?,
        None => bcf::Writer::from_stdout(&header, false, false)?,
    };

    let mut outobs = if let Some(f) = outobs {
        let mut writer = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(f)?;
        // write header for observations
        writer.write_record(
            [
                "chrom",
                "pos",
                "allele",
                "sample",
                "prob_mapping",
                "prob_mismapping",
                "prob_alt",
                "prob_ref",
                "prob_sample_alt",
                "evidence",
            ]
            .iter(),
        )?;
        Some(writer)
    } else {
        None
    };
    let mut record = inbcf.empty_record();
    let mut i = 0;
    loop {
        // read BCF
        if let Err(e) = inbcf.read(&mut record) {
            if e.is_eof() {
                return Ok(());
            } else {
                return Err(Box::new(e));
            }
        }
        i += 1;

        // translate to header of the writer
        outbcf.translate(&mut record);
        let mut pileups = pileups(
            &inbcf,
            &mut record,
            pair_model,
            &mut reference_buffer,
            omit_snvs,
            omit_indels,
            max_indel_len,
            exclusive_end,
        )?;

        if !pileups.is_empty() {
            // write observations
            if let Some(ref mut outobs) = outobs {
                let chrom = str::from_utf8(chrom(&inbcf, &record)).unwrap();
                for (i, pileup) in pileups.iter().enumerate() {
                    if let &Some(ref pileup) = pileup {
                        for obs in pileup.case_observations() {
                            outobs.serialize((chrom, record.pos(), i, "case", obs))?;
                        }
                        for obs in pileup.control_observations() {
                            outobs.serialize((chrom, record.pos(), i, "control", obs))?;
                        }
                    }
                }
                outobs.flush()?;
            }

            // write posterior probabilities
            let mut posterior_probs = Array::default((events.len(), pileups.len()));
            for (i, event) in events.iter().enumerate() {
                for (j, pileup) in pileups.iter_mut().enumerate() {
                    let p = if let &mut Some(ref mut pileup) = pileup {
                        // use joint probability instead of posterior since we do the
                        // normalization below, turning joint probabilities into posteriors.
                        Some(pileup.joint_prob(&event.af_case, &event.af_control))
                    } else {
                        // indicate missing value
                        None
                    };

                    posterior_probs[(i, j)] = p;
                }
            }

            for (j, pileup) in pileups.iter().enumerate() {
                if pileup.is_some() {
                    let total = LogProb::ln_sum_exp(
                        &posterior_probs
                            .column(j)
                            .iter()
                            .map(|v| v.unwrap())
                            .collect_vec(),
                    );
                    for i in 0..events.len() {
                        // normalize by total probability
                        let p = posterior_probs[(i, j)].unwrap();
                        posterior_probs[(i, j)] = Some(p - total);
                    }
                }
            }
            for (i, event) in events.iter().enumerate() {
                record.push_info_float(
                    event.tag_name("PROB").as_bytes(),
                    &phred_scale(posterior_probs.row(i).iter()),
                )?;
            }

            // write allele frequency estimates
            let mut case_afs = Vec::with_capacity(pileups.len());
            let mut control_afs = Vec::with_capacity(pileups.len());
            for pileup in &mut pileups {
                if let &mut Some(ref mut pileup) = pileup {
                    let (case_af, control_af) = pileup.map_allele_freqs();
                    case_afs.push(*case_af as f32);
                    control_afs.push(*control_af as f32);
                } else {
                    case_afs.push(f32::missing());
                    control_afs.push(f32::missing());
                }
            }
            record.push_info_float(b"CASE_AF", &case_afs)?;
            record.push_info_float(b"CONTROL_AF", &control_afs)?;
        }
        outbcf.write(&record)?;
        if i % 1000 == 0 {
            info!("{} records processed.", i);
        }
    }
}

use std::path::Path;
use std::error::Error;
use std::{f64, f32};
use std::str;

use itertools::Itertools;
use ndarray::prelude::*;
use csv;
use rust_htslib::bcf;
use bio::stats::{PHREDProb, LogProb};

use model::AlleleFreqs;
use model::priors;
use model::PairModel;
use model;
use ComplementEvent;
use Event;


const MISSING_VALUE: f32 = f32::NAN;


fn phred_scale<'a, I: IntoIterator<Item=&'a LogProb>>(probs: I) -> Vec<f32> {
    probs.into_iter().map(|&p| {
        if *p as f32 == MISSING_VALUE {
            MISSING_VALUE as f32
        } else {
            PHREDProb::from(p).abs() as f32
        }
    }).collect_vec()
}


pub struct PairEvent<A: AlleleFreqs, B: AlleleFreqs> {
    /// event name
    pub name: String,
    /// allele frequencies for case sample
    pub af_case: A,
    /// allele frequencies for control sample
    pub af_control: B
}


impl<A: AlleleFreqs, B: AlleleFreqs> Event for PairEvent<A, B> {
    fn name(&self) -> &str {
        &self.name
    }
}


fn pileups<'a, A, B, P>(inbcf: &bcf::Reader, record: &mut bcf::Record, joint_model: &'a mut PairModel<A, B, P>, omit_snvs: bool, omit_indels: bool) -> Result<Vec<Option<model::PairPileup<'a, A, B, P>>>, Box<Error>> where
    A: AlleleFreqs,
    B: AlleleFreqs,
    P: priors::PairModel<A, B>
{
    let chrom = chrom(&inbcf, &record);
    // TODO avoid cloning svtype
    let svtypes = record.info(b"SVTYPE").string().map(|values| {
        values.into_iter().map(|svtype| svtype.to_owned()).collect_vec()
    });

    let variants = if let Ok(svtypes) = svtypes {
        // obtain svlen if it is present or store error
        let svlens = record.info(b"SVLEN").integer().map(|values| values.to_owned());
        svtypes.iter().zip(try!(svlens)).map(|(svtype, svlen)| {
            if omit_indels {
                None
            } else if svtype == b"INS" {
                Some(model::Variant::Insertion(svlen.abs() as u32))
            } else if svtype == b"DEL" {
                Some(model::Variant::Deletion(svlen.abs() as u32))
            } else {
                None
            }
        }).collect_vec()
    } else {
        let alleles = record.alleles();
        let ref_allele = alleles[0];

        alleles.iter().skip(1).map(|alt_allele| {
            if alt_allele.len() == 1 && ref_allele.len() == 1 {
                if omit_snvs {
                    None
                } else {
                    Some(model::Variant::SNV(alt_allele[0]))
                }
            } else if alt_allele.len() == ref_allele.len() {
                // neither indel nor SNV
                None
            } else if omit_indels {
                None
            } else if alt_allele.len() < ref_allele.len() {
                Some(model::Variant::Deletion((ref_allele.len() - alt_allele.len()) as u32))
            } else {
                Some(model::Variant::Insertion((alt_allele.len() - ref_allele.len()) as u32))
            }
        }).collect_vec()
    };

    let mut pileups = Vec::with_capacity(variants.len());
    for variant in variants {
        pileups.push(if let Some(variant) = variant {
            Some(try!(joint_model.pileup(chrom, record.pos(), variant)))
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
/// * `inbcf` - path to BCF/VCF with preprocessed variant calls (`"-"` for STDIN).
/// * `outbcf` - path to BCF/VCF with results (`"-"` for STDOUT).
/// * `events` - events to call
/// * `complement_event` - optional complementary event to call (e.g. absent)
/// * `joint_model` - calling model to use
/// * `omit_snvs` - omit single nucleotide variants
/// * `omit_indels` - omit indels
/// * `outobs` - optional path where to store observations as JSON
///
/// # Returns
///
/// `Result` object with eventual error message.
pub fn call<A, B, P, M, R, W, X>(
    inbcf: &R,
    outbcf: &W,
    events: &[PairEvent<A, B>],
    complement_event: Option<&ComplementEvent>,
    pair_model: &mut PairModel<A, B, P>,
    omit_snvs: bool,
    omit_indels: bool,
    outobs: Option<&X>
) -> Result<(), Box<Error>> where
    A: AlleleFreqs,
    B: AlleleFreqs,
    P: priors::PairModel<A, B>,
    R: AsRef<Path>,
    W: AsRef<Path>,
    X: AsRef<Path>
{
    let inbcf = try!(bcf::Reader::new(inbcf));
    let mut header = bcf::Header::with_template(&inbcf.header);
    for event in events {
        header.push_record(
            event.header_entry("PROB", "PHRED-scaled probability for").as_bytes()
        );
    }
    if let Some(complement_event) = complement_event {
        header.push_record(complement_event.header_entry("PROB", "PHRED-scaled probability for").as_bytes());
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

    let mut outbcf = try!(bcf::Writer::new(outbcf, &header, false, false));
    let mut outobs = if let Some(f) = outobs {
        let mut writer = try!(csv::Writer::from_file(f)).delimiter(b'\t');
        // write header for observations
        try!(writer.write(["chrom", "pos", "allele", "sample", "prob_mapping", "prob_alt", "prob_ref", "prob_mismapped", "evidence"].iter()));
        Some(writer)
    } else { None };
    let mut record = bcf::Record::new();
    let mut i = 0;
    loop {
        if let Err(e) = inbcf.read(&mut record) {
            if e.is_eof() {
                return Ok(())
            } else {
                return Err(Box::new(e));
            }
        }
        i += 1;
        // translate to header of the writer
        outbcf.translate(&mut record);
        let pileups = try!(pileups(&inbcf, &mut record, pair_model, omit_snvs, omit_indels));

        if !pileups.is_empty() {
            if let Some(ref mut outobs) = outobs {
                let chrom = str::from_utf8(chrom(&inbcf, &record)).unwrap();
                for (i, pileup) in pileups.iter().enumerate() {
                    if let &Some(ref pileup) = pileup {
                        for obs in pileup.case_observations() {
                            try!(outobs.encode((chrom, record.pos(), i, "case", obs)));
                        }
                        for obs in pileup.control_observations() {
                            try!(outobs.encode((chrom, record.pos(), i, "control", obs)));
                        }
                    }
                }
                try!(outobs.flush());
            }

            let mut posterior_probs = Array::default((events.len(), pileups.len()));
            for (i, event) in events.iter().enumerate() {
                for (j, pileup) in pileups.iter().enumerate() {
                    let p = if let &Some(ref pileup) = pileup {
                        pileup.posterior_prob(&event.af_case, &event.af_control)
                    } else {
                        // indicate missing value
                        LogProb(MISSING_VALUE as f64)
                    };

                    posterior_probs[(i, j)] = p;
                }
                try!(record.push_info_float(
                    event.tag_name("PROB").as_bytes(),
                    &phred_scale(posterior_probs.row(i).iter())
                ));
            }
            if let Some(complement_event) = complement_event {
                let mut complement_probs = Vec::with_capacity(pileups.len());
                for (j, pileup) in pileups.iter().enumerate() {
                    let p = if pileup.is_some() {
                        let event_probs = posterior_probs.column(j).iter().cloned().collect_vec();
                        let total = LogProb::ln_sum_exp(&event_probs);
                        // total can slightly exceed 1 due to the numerical integration
                        if total > LogProb::ln_one() {
                            LogProb::ln_zero()
                        } else {
                            total.ln_one_minus_exp()
                        }
                    } else {
                        // indicate missing value
                        LogProb(MISSING_VALUE as f64)
                    };
                    complement_probs.push(p);
                }
                try!(record.push_info_float(
                    complement_event.tag_name("PROB").as_bytes(),
                    &phred_scale(complement_probs.iter())
                ));
            }
            let mut case_afs = Vec::with_capacity(pileups.len());
            let mut control_afs = Vec::with_capacity(pileups.len());
            for pileup in &pileups {
                if let &Some(ref pileup) = pileup {
                    let (case_af, control_af) = pileup.map_allele_freqs();
                    case_afs.push(*case_af as f32);
                    control_afs.push(*control_af as f32);
                } else {
                    case_afs.push(MISSING_VALUE);
                    control_afs.push(MISSING_VALUE);
                }
            }
            try!(record.push_info_float(b"CASE_AF", &case_afs));
            try!(record.push_info_float(b"CONTROL_AF", &control_afs));


        }
        try!(outbcf.write(&record));
        if i % 1000 == 0 {
            info!("{} records processed.", i);
        }
    }
}

fn chrom<'a>(inbcf: &'a bcf::Reader, record: &bcf::Record) -> &'a [u8] {
    inbcf.header.rid2name(record.rid().unwrap())
}

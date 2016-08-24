#![cfg_attr(feature="flame_it", feature(plugin))]
#![cfg_attr(feature="flame_it", plugin(flamer))]
// activate flame for the whole crate

extern crate bio;
extern crate rust_htslib;
#[macro_use]
extern crate log;
extern crate rgsl;
extern crate itertools;
#[macro_use]
extern crate approx;
extern crate rusty_machine;
extern crate ordered_float;
#[macro_use]
extern crate ndarray;
extern crate rustc_serialize;
#[macro_use]
extern crate quick_error;

#[cfg(feature="flame_it")]
extern crate flame;

pub mod model;
pub mod estimation;

pub use model::sample::Sample;
pub use model::likelihood;
pub use model::priors;
pub use model::sample::InsertSize;

use std::ascii::AsciiExt;


/// Event to call.
pub trait Event {
    fn name(&self) -> &str;

    fn tag_name(&self, prefix: &str) -> String {
        format!("{}_{}", prefix, self.name().to_ascii_uppercase())
    }

    fn header_entry(&self, prefix: &str, desc: &str) -> String {
        format!(
            "##INFO=<ID={tag_name},Number=A,Type=Float,\
            Description=\"{desc} {name} variant\">",
            name=self.name(),
            desc=desc,
            tag_name=&self.tag_name(prefix)
        )
    }
}


/// Complement of other given events (i.e. 1 - Pr(other events)).
pub struct ComplementEvent {
    /// event name
    pub name: String
}


impl Event for ComplementEvent {
    fn name(&self) -> &str {
        &self.name
    }
}


pub mod case_control {
    use std::path::Path;
    use std::error::Error;
    use std::{f64, f32};

    use itertools::Itertools;
    use ndarray::prelude::*;
    use rust_htslib::bcf;
    use bio::stats::{PHREDProb, LogProb};

    use model::AlleleFreqs;
    use model::priors;
    use model::JointModel;
    use model;
    use ComplementEvent;
    use Event;


    fn phred_scale<'a, I: IntoIterator<Item=&'a LogProb>>(probs: I) -> Vec<f32> {
        probs.into_iter().map(|&p| *PHREDProb::from(p) as f32).collect_vec()
    }


    pub struct CaseControlEvent<A: AlleleFreqs, B: AlleleFreqs> {
        /// event name
        pub name: String,
        /// allele frequencies for case sample
        pub af_case: A,
        /// allele frequencies for control sample
        pub af_control: B
    }


    impl<A: AlleleFreqs, B: AlleleFreqs> super::Event for CaseControlEvent<A, B> {
        fn name(&self) -> &str {
            &self.name
        }
    }


    fn pileups<A, B, P, M>(inbcf: &bcf::Reader, record: &mut bcf::Record, joint_model: &mut M, omit_snvs: bool, omit_indels: bool) -> Result<Vec<Option<model::Pileup<A, B, P>>>, Box<Error>> where
        A: AlleleFreqs,
        B: AlleleFreqs,
        P: priors::PairModel<A, B>,
        M: JointModel<A, B, P>
    {
        let chrom = inbcf.header.rid2name(record.rid().expect("reference id not found in header"));
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

            alleles.iter().map(|alt_allele| {
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
    /// * `events` - Events to call.
    /// * `joint_model` - Calling model to use.
    ///
    /// # Returns
    ///
    /// `Result` object with eventual error message.
    pub fn call<A, B, P, M, R, W>(
        inbcf: &R,
        outbcf: &W,
        events: &[CaseControlEvent<A, B>],
        complement_event: Option<&ComplementEvent>,
        joint_model: &mut M,
        omit_snvs: bool,
        omit_indels: bool
    ) -> Result<(), Box<Error>> where
        A: AlleleFreqs,
        B: AlleleFreqs,
        P: priors::PairModel<A, B>,
        M: JointModel<A, B, P>,
        R: AsRef<Path>,
        W: AsRef<Path>
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
            let pileups = try!(pileups(&inbcf, &mut record, joint_model, omit_snvs, omit_indels));

            if !pileups.is_empty() {
                let mut posterior_probs = Array::default((events.len(), pileups.len()));
                for (i, event) in events.iter().enumerate() {
                    for (j, pileup) in pileups.iter().enumerate() {
                        let p = if let &Some(ref pileup) = pileup {
                            pileup.posterior_prob(joint_model, &event.af_case, &event.af_control)
                        } else {
                            // indicate missing value
                            LogProb(f64::NAN)
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
                    for j in 0..pileups.len() {
                        let event_probs = posterior_probs.column(j).iter().cloned().collect_vec();
                        let total = LogProb::ln_sum_exp(&event_probs);
                        // total can slightly exceed 1 due to the numerical integration
                        let p = if total > LogProb::ln_one() {
                            LogProb::ln_zero()
                        } else {
                            total.ln_one_minus_exp()
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
                        let (case_af, control_af) = pileup.map_allele_freqs(joint_model);
                        case_afs.push(*case_af as f32);
                        control_afs.push(*control_af as f32);
                    } else {
                        case_afs.push(f32::NAN);
                        control_afs.push(f32::NAN);
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
}

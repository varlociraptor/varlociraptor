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

pub mod model;
pub mod estimation;

pub use model::sample::Sample;
pub use model::likelihood;
pub use model::priors;
pub use model::sample::InsertSize;


/// Event to call.
pub trait Event {
    fn name(&self) -> &str;
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
    use std::ascii::AsciiExt;
    use std::path::Path;
    use std::collections::HashMap;
    use itertools::Itertools;
    use rust_htslib::bcf;
    use bio::stats::logprobs;

    use model::priors::AlleleFreq;
    use model::priors;
    use model::JointModel;
    use model;
    use ComplementEvent;
    use Event;


    pub struct CaseControlEvent<A: AlleleFreq, B: AlleleFreq> {
        /// event name
        pub name: String,
        /// allele frequencies for case sample
        pub af_case: A,
        /// allele frequencies for control sample
        pub af_control: B
    }

    impl<A: AlleleFreq, B: AlleleFreq> super::Event for CaseControlEvent<A, B> {
        fn name(&self) -> &str {
            &self.name
        }
    }


    fn pileups<A, B, P, Q, M>(inbcf: &bcf::Reader, record: &mut bcf::Record, joint_model: &mut M) -> Result<Vec<model::Pileup<A, B, P, Q, M>>, String> where
        A: AlleleFreq,
        B: AlleleFreq,
        P: priors::Model<A>,
        Q: priors::Model<B>,
        M: JointModel<A, B, P, Q>
    {
        let svlen = record.info(b"SVLEN").integer().ok().map(|svlen| svlen[0]);
        let alleles = record.alleles();
        let mut pileups = Vec::with_capacity(alleles.len() - 1);
        let ref_allele = alleles[0];
        for alt_allele in &alleles {
            if alt_allele.len() == ref_allele.len() {
                // no indel
                continue;
            }

            let variant = if alt_allele == b"<DEL>" {
                if let Some(length) = svlen {
                    model::Variant::Deletion((length.abs()) as u32)
                } else {
                    return Err("Error reading SVLEN of <DEL> variant.".to_owned());
                }
            } else if alt_allele == b"<INS>" {
                if let Some(length) = svlen {
                    model::Variant::Insertion(length as u32)
                } else {
                    return Err("Error reading SVLEN of <INS> variant.".to_owned());
                }
            } else if alt_allele[0] == b'<' {
                // other special tag
                continue;
            } else if alt_allele.len() < ref_allele.len() {
                model::Variant::Deletion((ref_allele.len() - alt_allele.len()) as u32)
            } else {
                model::Variant::Insertion((alt_allele.len() - ref_allele.len()) as u32)
            };

            if let Some(rid) = record.rid() {
                let chrom = inbcf.header.rid2name(rid);
                // obtain pileup and calculate marginal probability
                let pileup = try!(joint_model.pileup(chrom, record.pos(), variant));
                pileups.push(pileup);
            } else {
                return Err("Error reading BCF/VCF record.".to_owned());
            }
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
    pub fn call<A, B, P, Q, M, R, W>(
        inbcf: &R,
        outbcf: &W,
        events: &[CaseControlEvent<A, B>],
        complement_event: Option<ComplementEvent>,
        joint_model: &mut M
    ) -> Result<(), String> where
        A: AlleleFreq,
        B: AlleleFreq,
        P: priors::Model<A>,
        Q: priors::Model<B>,
        M: JointModel<A, B, P, Q>,
        R: AsRef<Path>,
        W: AsRef<Path>
    {
        if let Ok(inbcf) = bcf::Reader::new(inbcf) {
            let mut header = bcf::Header::with_template(&inbcf.header);
            for event in events {
                header.push_record(
                    format!(
                        "##INFO=<ID=PROB_{name_upper},Number=A,Type=Float,\
                        Description=\"PHRED-scaled probability for {name} variant\">",
                        name=event.name(),
                        name_upper=&event.name().to_ascii_uppercase()
                    ).as_bytes()
                );
            }

            if let Ok(mut outbcf) = bcf::Writer::new(outbcf, &header, false, false) {
                let mut record = bcf::Record::new();
                let mut i = 0;
                loop {
                    // Read BCF/VCF record.
                    match inbcf.read(&mut record) {
                        Err(bcf::ReadError::NoMoreRecord) => break,
                        Err(_) => return Err("Error reading BCF/VCF.".to_owned()),
                        _ => ()
                    }
                    i += 1;
                    // translate to header of the writer
                    outbcf.translate(&mut record);
                    let pileups = try!(pileups(&inbcf, &mut record, joint_model));

                    let mut posterior_probs = Vec::with_capacity(record.alleles().len() - 1);
                    for event in events {
                        posterior_probs.clear();
                        for pileup in &pileups {
                            let p = pileup.posterior_prob(&event.af_case, &event.af_control);
                            posterior_probs.push(logprobs::log_to_phred(p) as f32);
                        }
                        if !posterior_probs.is_empty() {
                            if record.push_info_float(
                                format!("PROB_{}", event.name().to_ascii_uppercase()).as_bytes(),
                                &posterior_probs
                            ).is_err() {
                                return Err("Error writing INFO tag in BCF.".to_owned());
                            }
                        }
                    }
                    // TODO handle complement

                    if outbcf.write(&record).is_err() {
                        return Err("Error writing BCF record.".to_owned());
                    }
                    if i % 1000 == 0 {
                        info!("{} records processed.", i);
                    }
                }
                Ok(())
            } else {
                Err("Error writing BCF record.".to_owned())
            }
        } else {
            Err("Error reading BCF/VCF".to_owned())
        }
    }
}

extern crate bio;
extern crate rust_htslib;
#[macro_use]
extern crate log;
extern crate rgsl;
extern crate itertools;

pub mod model;

use std::ops::Range;
use rust_htslib::bcf;
use bio::stats::logprobs;

use model::priors;


/// Event to call.
pub struct Event {
    /// BCF/VCF to store results
    pub writer: bcf::Writer,
    /// continuous allele frequency range (case sample)
    pub af_case: Range<f64>,
    /// discrete allele frequencies (control sample)
    pub af_control: Vec<f64>
}


/// Call variants with the given model.
///
/// # Arguments
///
/// * `bcf` - BCF/VCF reader with preprocessed variant calls.
/// * `events` - Events to call.
/// * `joint_model` - Calling model to use.
///
/// # Returns
///
/// `Result` object with eventual error message.
pub fn call<P: priors::ContinuousModel, Q: priors::DiscreteModel>(
    bcf: &mut bcf::Reader,
    events: &mut [Event],
    joint_model: &mut model::JointModel<P, Q>
) -> Result<(), String> {

    let mut record = bcf::Record::new();
    loop {
        // Read BCF/VCF record.
        match bcf.read(&mut record) {
            Err(bcf::ReadError::NoMoreRecord) => break,
            Err(_) => return Err("Error reading BCF/VCF.".to_owned()),
            _ => ()
        }

        // Iterate over alleles.
        let ref_allele = record.alleles()[0].to_owned();
        for i in 1..record.allele_count() as usize {
            let alt_allele = record.alleles()[i].to_owned();
            if alt_allele.len() == ref_allele.len() {
                // no indel
                continue;
            }
            let is_del = alt_allele.len() < ref_allele.len();
            if let Some(rid) = record.rid() {
                let chrom = bcf.header.rid2name(rid);
                // obtain pileup and calculate marginal probability
                let pileup = try!(joint_model.pileup(chrom, record.pos(), alt_allele.len() as u32, is_del));
                for event in events.iter_mut() {
                    // calculate posterior probability for event
                    let posterior_prob = try!(pileup.posterior_prob(&event.af_case, &event.af_control));
                    // TODO properly handle multiple alleles!
                    record.set_qual(logprobs::log_to_phred(posterior_prob) as f32);
                    if let Err(_) = event.writer.write(&record) {
                        return Err("Error writing BCF/VCF record.".to_owned());
                    }
                }
            } else {
                return Err("Error reading BCF/VCF.".to_owned());
            }
        }
    }
    Ok(())
}

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

use model::sample::Sample;


pub struct CallClass {
    writer: bcf::Writer,
    af_case: Range<f64>,
    af_control: Vec<f64>
}


pub fn call<A: Sample, B: Sample>(
    bcf: &mut bcf::Reader,
    classes: &mut [CallClass],
    joint_model: &mut model::JointModel<A, B>
) -> Result<(), String> {

    for record in bcf.records() {
        if let Ok(mut record) = record {
            let ref_allele = record.alleles()[0].to_owned();
            for i in 1..record.allele_count() as usize {
                let alt_allele = record.alleles()[i].to_owned();
                if alt_allele.len() == ref_allele.len() {
                    continue;
                }
                let is_del = alt_allele.len() < ref_allele.len();
                if let Some(rid) = record.rid() {
                    let chrom = bcf.header.rid2name(rid);
                    let pileup = try!(joint_model.pileup(chrom, record.pos(), alt_allele.len() as u32, is_del));
                    for class in classes.iter_mut() {
                        let posterior_prob = try!(pileup.posterior_prob(&class.af_case, &class.af_control));
                        record.set_qual(logprobs::log_to_phred(posterior_prob) as f32);
                        if let Err(_) = class.writer.write(&record) {
                            return Err("Error writing BCF/VCF record.".to_owned());
                        }
                    }
                } else {
                    return Err("Error reading BCF/VCF.".to_owned());
                }
            }
        } else {
            return Err("Error reading VCF".to_owned());
        }
    }
    Ok(())
}

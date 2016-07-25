extern crate bio;
extern crate rust_htslib;
#[macro_use]
extern crate log;
extern crate rgsl;
extern crate itertools;

pub mod model;

use std::ops::Range;
use rust_htslib::bam;
use rust_htslib::bcf;
use bio::stats::logprobs;

use model::observations;
use model::InsertSize;
use model::priors;


pub struct CallClass {
    writer: bcf::Writer,
    af_case: Range<f64>,
    af_control: Vec<f64>
}


pub fn call<P: priors::Model, Q: priors::Model>(
    case_bam: bam::IndexedReader,
    control_bam: bam::IndexedReader,
    bcf: &mut bcf::Reader,
    classes: &mut [CallClass],
    joint_model: model::JointModel<P, Q>,
    case_insert_size: InsertSize,
    control_insert_size: InsertSize,
    pileup_window: u32
) -> Result<(), String> {
    let mut case_processor = observations::BAMProcessor::new(case_bam, case_insert_size, pileup_window);
    let mut control_processor = observations::BAMProcessor::new(control_bam, control_insert_size, pileup_window);

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
                    let case_pileup = try!(case_processor.extract_observations(chrom, record.pos(), alt_allele.len() as u32, is_del));
                    let control_pileup = try!(control_processor.extract_observations(chrom, record.pos(), alt_allele.len() as u32, is_del));
                    let marginal_prob = try!(joint_model.marginal_prob(&case_pileup, &control_pileup));
                    for class in classes.iter_mut() {
                        let posterior_prob = try!(joint_model.posterior_prob(&case_pileup, &control_pileup, &class.af_case, &class.af_control, marginal_prob));
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

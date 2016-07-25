use rust_htslib::bam;
use rust_htslib::bcf;

pub mod likelihood;
pub mod observations;


/*
pub fn call_tumor_normal(tumor_bam: bam::IndexedReader, normal_bam: bam::IndexedReader, in_bcf: bcf::Reader, out_bcf: bcf::Writer, pileup_window: u32) {
    let bam_processor = observations::BAMProcessor::new(tumor_bam, pileup_window);
    let bam_processor = observations::BAMProcessor::new(tumor_bam, pileup_window);
    for record in in_bcf.records() {

    }
}*/

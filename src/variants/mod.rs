use std::ops::Deref;

use anyhow::Result;
use bio::stats::LogProb;
use rust_htslib::bam;

use crate::utils::{Overlap, GenomicLocus};
use crate::model::sample;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::model::evidence::{observation::ObservationBuilder, Observation};

pub mod mnv;
pub mod snv;

#[derive(Debug, Getters, new)]
#[get = "pub"]
pub struct AlleleProb {
    ref_allele: LogProb,
    alt_allele: LogProb,
}

pub trait Variant<'a> {
    type Evidence: Evidence<'a>;
    type Loci: Loci;

    /// Determine whether record overlaps the variant.
    fn overlap(&self, evidence: &Self::Evidence) -> Overlap;

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci;

    /// Calculate probability for alt and reference allele.
    fn prob_alleles(&self, evidence: &Self::Evidence) -> Result<Option<AlleleProb>>;

    /// Calculate probability to sample a record length like the given one from the alt allele.
    fn prob_sample_alt(&self, evidence: &Self::Evidence, alignment_properties: &AlignmentProperties) -> LogProb;
}

pub trait Observable {
    fn extract_observations(&self, buffer: &mut sample::RecordBuffer) -> Result<Vec<Observation>>;
}

impl<'a, V> Observable for V where V: Variant<'a, Evidence=SingleEndEvidence<'a>, Loci=SingleLocus> {
    fn extract_observations(&self, buffer: &mut sample::RecordBuffer) -> Result<Vec<Observation>> {
        let locus = self.loci();
        buffer.fetch_surrounding(locus.chrom(), locus.pos(), locus.pos() + locus.len())?;

        let candidates = buffer.iter().filter_map(|record| {
            let evidence = SingleEndEvidence::new(record);
            if record.pos() as u32 > locus.pos() || !self.overlap(&evidence).is_enclosing() {
                None
            } else {
                Some(evidence)
            }
        }).collect();

        
    }
}


// pub trait SmallVariant {
    
// }

// pub trait StructuralVariant: Variant {
//     /// Calculate probability for alt and reference allele.
//     fn prob_alleles(
//         &self,
//         left: &mut bam::Record,
//         right: &mut bam::Record,
//     ) -> Result<Option<AlleleProb>>;

//     /// Calculate probability to sample a record length like the given one from the alt allele.
//     fn prob_sample_alt(&self, left: &mut bam::Record, right: &mut bam::Record, alignment_properties: &AlignmentProperties) -> LogProb;
// }

pub trait Evidence<'a> {}

#[derive(new)]
pub struct SingleEndEvidence<'a> {
    inner: &'a bam::Record
}

impl<'a> Deref for SingleEndEvidence<'a> {
    type Target = bam::Record;

    fn deref(&self) -> &bam::Record {
        self.inner
    }
}

pub enum PairedEndEvidence<'a> {
    SingleEnd(&'a bam::Record),
    PairedEnd{ left: &'a bam::Record, right: &'a bam::Record }
}

impl<'a> Evidence<'a> for SingleEndEvidence<'a> {}

impl<'a> Evidence<'a> for PairedEndEvidence<'a> {}


pub trait Loci {}

#[derive(new, Debug, Getters, CopyGetters, Derefable)]
pub struct SingleLocus {
    #[getset(get = "pub")]
    #[deref]
    locus: GenomicLocus,
    #[getset(get_copy = "pub")]
    len: u32
}

impl Loci for SingleLocus {}

#[derive(new, Debug)]
pub struct MultiLocus {
    loci: Vec<SingleLocus>,
}

impl Loci for MultiLocus {}
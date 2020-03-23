use anyhow::Result;
use bio::stats::LogProb;
use rust_htslib::bam;

use crate::utils::Overlap;

pub mod mnv;
pub mod snv;

#[derive(Debug, Getters, new)]
#[get = "pub"]
pub struct AlleleProb {
    ref_allele: LogProb,
    alt_allele: LogProb,
}

pub trait Variant {
    /// Determine whether read overlaps the variant.
    fn overlap(&self, read: &mut bam::Record) -> Overlap;
}

pub trait ReadEvidence {
    /// Calculate probability for alt and reference allele.
    fn prob_alleles(&self, read: &mut bam::Record) -> Result<Option<AlleleProb>>;

    /// Calculate probability to sample a read length like the given one from the alt allele.
    fn prob_sample_alt(&self, read: &mut bam::Record) -> LogProb;
}

pub trait ReadPairEvidence {
    /// Calculate probability for alt and reference allele.
    fn prob_alleles(
        &self,
        left: &mut bam::Record,
        right: &mut bam::Record,
    ) -> Result<Option<AlleleProb>>;

    /// Calculate probability to sample a read length like the given one from the alt allele.
    fn prob_sample_alt(&self, left: &mut bam::Record, right: &mut bam::Record) -> LogProb;
}


pub trait EvidenceFetcher {
    fn fetch(buffer: &mut bam::buffer::RecordBuffer) -> Result<()>;
}

pub struct SingleReadEvidenceFetcher<V: Variant> {
    reads: Vec<bam::Record>,
    variant: V
}

impl EvidenceFetcher for SingleReadEvidenceFetcher {
    fn fetch(buffer: &mut bam::buffer::RecordBuffer) -> Result<()> {
        buffer.fetch(chrom: &[u8], start: u32, end: u32)
    }
}
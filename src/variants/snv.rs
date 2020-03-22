use rust_htslib::bam;
use bio::stats::LogProb;
use anyhow::Result;

use crate::variants::{Variant, AlleleProb, ReadEvidence};
use crate::utils::{Overlap, GenomicLocus};
use crate::model::evidence::reads::prob_read_base;

pub struct SNV {
    locus: GenomicLocus,
    ref_base: u8,
    alt_base: u8,
}

impl Variant for SNV {
    fn overlap(&self, read: &mut bam::Record) -> Overlap {
        let read_start = read.pos() as u32;
        let read_end = read.cigar_cached().unwrap().end_pos() as u32;
        if read_start <= self.locus.pos() && read_end > self.locus.pos() {
            Overlap::Enclosing(1)
        } else {
            Overlap::None
        }
    }
}

impl ReadEvidence for SNV {
    fn prob_alleles(&self, read: &mut bam::Record) -> Result<Option<AlleleProb>> {
        if let Some(qpos) = read.cigar_cached().unwrap().read_pos(self.locus.pos(), false, false)? {
            let read_base = read.seq()[qpos as usize];
            let base_qual = read.qual()[qpos as usize];
            let prob_alt = prob_read_base(read_base, self.alt_base, base_qual);

            // METHOD: instead of considering the actual REF base, we assume that REF is whatever
            // base the read has at this position (if not the ALT base). This way, we avoid biased
            // allele frequencies at sites with multiple alternative alleles.
            // Note that this is an approximation. The real solution would be to have multiple allele
            // frequency variables in the likelihood function, but that would be computationally
            // more demanding (leading to a combinatorial explosion).
            // However, the approximation is pretty accurate, because it will only matter for true
            // multiallelic cases. Sequencing errors won't have a severe effect on the allele frequencies
            // because they are too rare.
            let non_alt_base = if read_base != self.alt_base {
                read_base
            } else {
                self.ref_base
            };

            let prob_ref = prob_read_base(read_base, non_alt_base, base_qual);
            Ok(Some(AlleleProb::new(prob_ref, prob_alt)))
        } else {
            // a read that spans an SNV might have the respective position deleted (Cigar op 'D')
            // or reference skipped (Cigar op 'N'), and the library should not choke on those reads
            // but instead needs to know NOT to add those reads (as observations) further up
            Ok(None)
        }
    }

    fn prob_sample_alt(&self, read: &mut bam::Record) -> LogProb {
        LogProb::ln_one()
    }
}
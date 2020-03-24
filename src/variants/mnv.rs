use anyhow::Result;
use bio::stats::LogProb;

use crate::model::evidence::reads::prob_read_base;
use crate::utils::{Overlap, GenomicLocus};
use crate::variants::{AlleleProb, Variant, SingleEndEvidence, SingleLocus};
use crate::estimation::alignment_properties::AlignmentProperties;

pub struct MNV {
    locus: SingleLocus,
    ref_bases: Vec<u8>,
    alt_bases: Vec<u8>,
}

impl MNV {
    pub fn new(locus: GenomicLocus, ref_bases: Vec<u8>, alt_bases: Vec<u8>) -> Self {
        MNV {
            locus: SingleLocus::new(locus, alt_bases.len() as u32),
            ref_bases,
            alt_bases,
        }
    }

    pub fn locus(&self) -> &SingleLocus {
        self.loci()
    }
}

impl<'a> Variant<'a> for MNV {
    type Evidence = SingleEndEvidence<'a>;
    type Loci = SingleLocus;

    fn overlap(&self, read: &SingleEndEvidence) -> Overlap {
        let read_start = read.pos() as u32;
        let read_end = read.cigar_cached().unwrap().end_pos() as u32;
        if read_start <= self.locus.pos() && read_end > self.locus.pos() + self.locus().len() {
            Overlap::Enclosing(self.locus().len())
        } else {
            Overlap::None
        }
    }

    fn loci(&self) -> &SingleLocus {
        &self.locus
    }

    fn prob_alleles(&self, read: &SingleEndEvidence) -> Result<Option<AlleleProb>> {
        let mut prob_ref = LogProb::ln_one();
        let mut prob_alt = LogProb::ln_one();
        for ((alt_base, ref_base), pos) in self
            .alt_bases
            .iter()
            .zip(self.ref_bases.iter())
            .zip(self.locus.pos()..self.locus.pos() + self.locus().len())
        {
            if let Some(qpos) = read.cigar_cached().unwrap().read_pos(pos, false, false)? {
                let read_base = read.seq()[qpos as usize];
                let base_qual = read.qual()[qpos as usize];

                // METHOD: instead of considering the actual REF base, we assume that REF is whatever
                // base the read has at this position (if not the ALT base). This way, we avoid biased
                // allele frequencies at sites with multiple alternative alleles.
                // Note that this is an approximation. The real solution would be to have multiple allele
                // frequency variables in the likelihood function, but that would be computationally
                // more demanding (leading to a combinatorial explosion).
                // However, the approximation is pretty accurate, because it will only matter for true
                // multiallelic cases. Sequencing errors won't have a severe effect on the allele frequencies
                // because they are too rare.
                let non_alt_base = if read_base != *alt_base {
                    read_base
                } else {
                    *ref_base
                };

                prob_alt += prob_read_base(read_base, *alt_base, base_qual);
                prob_ref += prob_read_base(read_base, non_alt_base, base_qual);
            } else {
                // a read that spans an SNV might have the respective position deleted (Cigar op 'D')
                // or reference skipped (Cigar op 'N'), and the library should not choke on those reads
                // but instead needs to know NOT to add those reads (as observations) further up
                return Ok(None);
            }
        }
        Ok(Some(AlleleProb::new(prob_ref, prob_alt)))
    }

    fn prob_sample_alt(&self, _: &SingleEndEvidence, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

use anyhow::Result;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::model::evidence::reads::prob_read_base;
use crate::variants::{realignable, AlleleProb, Overlap, SingleEndEvidence, SingleLocus, Variant};

pub struct SNV {
    locus: SingleLocus,
    ref_base: u8,
    alt_base: u8,
}

impl SNV {
    pub fn new(locus: genome::Locus, ref_base: u8, alt_base: u8) -> Self {
        SNV {
            locus: SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + 1,
            )),
            ref_base,
            alt_base,
        }
    }
}

impl<'a> Variant<'a> for SNV {
    type Evidence = SingleEndEvidence<'a>;
    type Loci = SingleLocus;

    fn is_valid_evidence(&self, evidence: &SingleEndEvidence) -> Option<Vec<usize>> {
        if let Overlap::Enclosing = self.locus.overlap(evidence, false) {
            Some(vec![0])
        } else {
            None
        }
    }

    fn loci(&self) -> &SingleLocus {
        &self.locus
    }

    fn prob_alleles(
        &self,
        read: &SingleEndEvidence,
        _: &AlignmentProperties,
    ) -> Result<Option<AlleleProb>> {
        if let Some(qpos) = read
            .cigar_cached()
            .unwrap()
            // TODO expect u64 in read_pos
            .read_pos(self.locus.range().start as u32, false, false)?
        {
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

    fn prob_sample_alt(&self, _: &SingleEndEvidence, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

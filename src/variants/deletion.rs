use std::cmp;

use anyhow::Result;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval};
use rgsl::randist::gaussian::ugaussian_P;
use bio::stats::pairhmm;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::utils::NUMERICAL_EPSILON;
use crate::variants::{AlleleProb, MultiLocus, PairedEndEvidence, SingleLocus, Variant, FragmentEnclosable};
use crate::variants::realignable::{pairhmm::RefBaseEmission, pairhmm::default_ref_base_emission};

pub struct Deletion {
    locus: genome::Interval,
    fetch_loci: MultiLocus,
}

impl Deletion {
    pub fn new(locus: genome::Interval) -> Self {
        let start = locus.range().start;
        let end = locus.range().end;
        let centerpoint = start + ((end - start) as f64 / 2.0).round() as u64;
        let contig = locus.contig().to_owned();

        let fetch_loci = MultiLocus::new(vec![
            SingleLocus(genome::Interval::new(contig.clone(), start..start + 1)),
            SingleLocus(genome::Interval::new(
                contig.clone(),
                centerpoint..centerpoint + 1,
            )),
            SingleLocus(genome::Interval::new(contig, end - 1..end)),
        ]);

        Deletion {
            locus: locus,
            fetch_loci,
        }
    }
}

impl<'a> FragmentEnclosable<'a> for Deletion {
    fn len(&self) -> u64 {
        self.locus.range().end - self.locus.range().start
    }
}

impl<'a> Variant<'a> for Deletion {
    type Evidence = PairedEndEvidence<'a>;
    type Loci = MultiLocus;

    fn is_valid_evidence(&self, evidence: &Self::Evidence) -> Option<Vec<usize>> {
        let mut locus_idx = Vec::new();
        for (i, locus) in self.loci().iter().enumerate() {
            if match evidence {
                PairedEndEvidence::SingleEnd(read) => { !locus.overlap(read, true).is_none() }
                PairedEndEvidence::PairedEnd { left, right } => {
                    !locus.overlap(left, true).is_none() || !locus.overlap(right, true).is_none()
                }
            } {
                locus_idx.push(i);
            }
        }
        if locus_idx.is_empty() {
            None
        } else {
            Some(locus_idx)
        }
    }

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci {
        &self.fetch_loci
    }

    /// Calculate probability for alt and reference allele.
    fn prob_alleles(&self, evidence: &Self::Evidence) -> Result<Option<AlleleProb>> {

    }

    fn prob_sample_alt(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        match evidence {
            PairedEndEvidence::PairedEnd { left, right } => {
                self.prob_sample_alt_paired(left.seq().len() as u64, right.seq().len() as u64, alignment_properties)
            },
            PairedEndEvidence::SingleEnd( read ) => {
                self.prob_sample_alt_single(read.seq().len() as u64, alignment_properties)
            }
        }
    }
}


/// Emission parameters for PairHMM over deletion allele.
#[derive(Debug)]
pub struct DeletionEmissionParams<'a> {
    ref_seq: &'a [u8],
    ref_offset: usize,
    ref_end: usize,
    del_start: usize,
    del_len: usize,
    read_emission: &'a ReadEmission<'a>,
}

impl<'a> RefBaseEmission for DeletionEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;
        if i_ <= self.del_start {
            self.ref_seq[i_]
        } else {
            self.ref_seq[i_ + self.del_len]
        }
    }

    default_ref_base_emission!();
}

impl<'a> pairhmm::EmissionParameters for DeletionEmissionParams<'a> {
    default_emission!();

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }
}
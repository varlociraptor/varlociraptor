use std::cmp;
use std::sync::Arc;

use anyhow::Result;
use bio::stats::pairhmm::EmissionParameters;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval};
use rust_htslib::bam;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::insert_size::estimate_insert_size;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::sampling_bias::{FragmentSamplingBias, ReadSamplingBias, SamplingBias};
use crate::variants::{AlleleProb, MultiLocus, PairedEndEvidence, SingleLocus, Variant};
use crate::{default_emission, default_ref_base_emission};

pub struct Deletion {
    locus: genome::Interval,
    fetch_loci: MultiLocus,
    realigner: Realigner,
}

impl Deletion {
    pub fn new(locus: genome::Interval, realigner: Realigner) -> Self {
        let start = locus.range().start;
        let end = locus.range().end;
        let len = end - start;
        let centerpoint = start + (len as f64 / 2.0).round() as u64;
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
            realigner,
        }
    }

    pub fn centerpoint(&self) -> u64 {
        self.fetch_loci[1].range().start
    }

    pub fn prob_alleles_isize(
        &self,
        left_record: &bam::Record,
        right_record: &bam::Record,
        alignment_properties: &AlignmentProperties,
    ) -> Result<AlleleProb> {
        let insert_size = estimate_insert_size(left_record, right_record)?;

        let p_ref = self.isize_pmf(insert_size, 0.0, alignment_properties);
        let p_alt = self.isize_pmf(insert_size, self.len() as f64, alignment_properties);

        if (p_ref == LogProb::ln_zero()
            && !self.is_within_sd(insert_size, self.len() as f64, alignment_properties))
            || (p_alt == LogProb::ln_zero()
                && !self.is_within_sd(insert_size, 0.0, alignment_properties))
        {
            // METHOD: We cannot consider insert size as a reliable estimate here, because it is
            // outside of the numerical resolution for one of the alleles, and not within a
            // standard deviation away from the mean for the other allele.
            Ok(AlleleProb::new(LogProb::ln_one(), LogProb::ln_one()))
        } else {
            Ok(AlleleProb::new(p_ref, p_alt))
        }
    }
}

impl<'a> SamplingBias<'a> for Deletion {
    fn len(&self) -> u64 {
        self.locus.range().end - self.locus.range().start
    }
}

impl<'a> FragmentSamplingBias<'a> for Deletion {}
impl<'a> ReadSamplingBias<'a> for Deletion {}

impl<'a: 'c, 'b: 'c, 'c> Realignable<'a, 'b, 'c> for Deletion {
    type EmissionParams = DeletionEmissionParams<'c>;

    fn alt_emission_params(
        &'b self,
        read_emission_params: &'a ReadEmission,
        ref_seq: Arc<Vec<u8>>,
        ref_window: usize,
    ) -> DeletionEmissionParams<'c> {
        let start = self.locus.range().start as usize;
        let end = self.locus.range().end as usize;
        DeletionEmissionParams::<'c> {
            del_start: start,
            del_len: end - start,
            ref_offset: start.saturating_sub(ref_window),
            ref_end: cmp::min(start + ref_window, ref_seq.len() - self.len() as usize),
            ref_seq,
            read_emission: read_emission_params,
        }
    }
}

impl<'a> Variant<'a> for Deletion {
    type Evidence = PairedEndEvidence<'a>;
    type Loci = MultiLocus;

    fn is_valid_evidence(&self, evidence: &Self::Evidence) -> Option<Vec<usize>> {
        let mut locus_idx = Vec::new();
        for (i, locus) in self.loci().iter().enumerate() {
            if match evidence {
                PairedEndEvidence::SingleEnd(read) => !locus.overlap(read, true).is_none(),
                PairedEndEvidence::PairedEnd { left, right } => {
                    // METHOD: valid if one read overlaps and the fragment encloses the centerpoint.
                    let right_cigar = right.cigar_cached().unwrap();
                    (!locus.overlap(left, true).is_none() || !locus.overlap(right, true).is_none())
                        && (left.pos() as u64) < self.centerpoint()
                        && right_cigar.end_pos() as u64 > self.centerpoint()
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

    fn prob_alleles(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> Result<Option<AlleleProb>> {
        match evidence {
            PairedEndEvidence::SingleEnd(record) => Ok(Some(self.realigner.prob_alleles(
                record,
                &self.locus,
                self,
            )?)),
            PairedEndEvidence::PairedEnd { left, right } => {
                // METHOD: Extract insert size information for fragments (e.g. read pairs) spanning an indel of interest
                // Here we calculate the product of insert size based and alignment based probabilities.
                // This has the benefit that the calculation automatically checks for consistence between
                // insert size and overlapping alignmnments.
                // This sports the following desirable behavior:
                //
                // * If there is no clear evidence from either the insert size or the alignment, the factors
                //   simply scale because the probabilities of the corresponding type of evidence will be equal.
                // * If reads and fragments agree, 1 stays 1 and 0 stays 0.
                // * If reads and fragments disagree (the promising part!), the other type of evidence will
                //   scale potential false positive probabilities down.
                // * Since there is only one observation per fragment, there is no double counting when
                //   estimating allele frequencies. Before, we had one observation for an overlapping read
                //   and potentially another observation for the corresponding fragment.
                let prob_left = self.realigner.prob_alleles(left, &self.locus, self)?;
                let prob_right = self.realigner.prob_alleles(right, &self.locus, self)?;
                let prob_isize = self.prob_alleles_isize(left, right, alignment_properties)?;

                Ok(Some(AlleleProb::new(
                    prob_left.ref_allele() + prob_right.ref_allele() + prob_isize.ref_allele(),
                    prob_left.alt_allele() + prob_right.alt_allele() + prob_isize.alt_allele(),
                )))
            }
        }
    }

    fn prob_sample_alt(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        match evidence {
            PairedEndEvidence::PairedEnd { left, right } => self.prob_sample_alt_fragment(
                left.seq().len() as u64,
                right.seq().len() as u64,
                alignment_properties,
            ),
            PairedEndEvidence::SingleEnd(read) => {
                self.prob_sample_alt_read(read.seq().len() as u64, alignment_properties)
            }
        }
    }
}

/// Emission parameters for PairHMM over deletion allele.
#[derive(Debug)]
pub struct DeletionEmissionParams<'a> {
    ref_seq: Arc<Vec<u8>>,
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

impl<'a> EmissionParameters for DeletionEmissionParams<'a> {
    default_emission!();

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }
}

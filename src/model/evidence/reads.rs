// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::error::Error;
use std::fmt::Debug;

use bio::stats::{LogProb, PHREDProb, Prob};
use bio_types::strand::Strand;
use rust_htslib::bam;
use rust_htslib::bam::record::CigarStringView;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::model::Variant;
use bio::pattern_matching::myers::Myers;
use bio::stats::pairhmm;

/// Convert MAPQ (from read mapper) to LogProb for the event that the read maps
/// correctly and the event that it maps incorrectly.
fn prob_mapping_mismapping(record: &bam::Record) -> (LogProb, LogProb) {
    let prob_mismapping = LogProb::from(PHREDProb(record.mapq() as f64));
    let prob_mapping = prob_mismapping.ln_one_minus_exp();
    (prob_mapping, prob_mismapping)
}

pub trait AbstractReadEvidence {
    /// Calculate probability for reference and alternative allele.
    fn prob(
        &mut self,
        record: &bam::Record,
        cigar: &CigarStringView,
        start: u32,
        variant: &Variant,
        ref_seq: &[u8],
    ) -> Result<Option<(LogProb, LogProb)>, Box<dyn Error>>;

    /// Calculate mapping and mismapping probability of given record.
    fn prob_mapping_mismapping(&self, record: &bam::Record) -> (LogProb, LogProb) {
        prob_mapping_mismapping(record)
    }

    fn prob_sample_alt(
        &self,
        read_len: u32,
        variant: &Variant,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb;

    fn strand(&self, record: &bam::Record) -> Strand {
        if record.flags() & 0x10 != 0 {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }
}

#[derive(Debug, Clone)]
pub struct NoneEvidence;

impl NoneEvidence {
    pub fn new() -> Self {
        NoneEvidence
    }
}

impl AbstractReadEvidence for NoneEvidence {
    fn prob(
        &mut self,
        record: &bam::Record,
        cigar: &CigarStringView,
        start: u32,
        variant: &Variant,
        ref_seq: &[u8],
    ) -> Result<Option<(LogProb, LogProb)>, Box<dyn Error>> {
        // TODO: Should we make this check against potential indel alt alleles, as well? Would need to collect respective observations / reads, then.
        if let &Variant::None = variant {
            if let Some(qpos) = cigar.read_pos(start, false, false)? {
                let read_base = record.seq()[qpos as usize];
                let base_qual = record.qual()[qpos as usize];
                let miscall = prob_read_base_miscall(base_qual);
                // here, prob_alt is the probability of any alternative allele / nucleotide, NOT of a particular alternative allele
                if read_base.to_ascii_uppercase() == ref_seq[start as usize].to_ascii_uppercase() {
                    Ok(Some((miscall.ln_one_minus_exp(), miscall)))
                } else {
                    Ok(Some((miscall, miscall.ln_one_minus_exp())))
                }
            } else {
                // a read that spans a potential Ref site might have the respective position deleted (Cigar op 'D')
                // or reference skipped (Cigar op 'N'), and the library should not choke on those reads
                // but instead needs to know NOT to add those reads (as observations) further up
                Ok(None)
            }
        } else {
            panic!("bug: unsupported variant");
        }
    }

    fn prob_sample_alt(&self, _: u32, _: &Variant, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

#[derive(Debug, Clone)]
pub struct SNVEvidence;

impl SNVEvidence {
    pub fn new() -> Self {
        SNVEvidence
    }
}

impl AbstractReadEvidence for SNVEvidence {
    fn prob(
        &mut self,
        record: &bam::Record,
        cigar: &CigarStringView,
        start: u32,
        variant: &Variant,
        ref_seq: &[u8],
    ) -> Result<Option<(LogProb, LogProb)>, Box<dyn Error>> {
        if let &Variant::SNV(base) = variant {
            if let Some(qpos) = cigar.read_pos(start, false, false)? {
                let read_base = record.seq()[qpos as usize];
                let base_qual = record.qual()[qpos as usize];
                let prob_alt = prob_read_base(read_base, base, base_qual);
                let prob_ref = prob_read_base(read_base, ref_seq[start as usize], base_qual);
                Ok(Some((prob_ref, prob_alt)))
            } else {
                // a read that spans an SNV might have the respective position deleted (Cigar op 'D')
                // or reference skipped (Cigar op 'N'), and the library should not choke on those reads
                // but instead needs to know NOT to add those reads (as observations) further up
                Ok(None)
            }
        } else {
            panic!("bug: unsupported variant");
        }
    }

    fn prob_sample_alt(&self, _: u32, _: &Variant, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

#[derive(Debug, Clone)]
pub struct MNVEvidence;

impl MNVEvidence {
    pub fn new() -> Self {
        MNVEvidence
    }
}

impl AbstractReadEvidence for MNVEvidence {
    fn prob(
        &mut self,
        record: &bam::Record,
        cigar: &CigarStringView,
        start: u32,
        variant: &Variant,
        ref_seq: &[u8],
    ) -> Result<Option<(LogProb, LogProb)>, Box<dyn Error>> {
        if let &Variant::MNV(ref bases) = variant {
            let mut prob_ref = LogProb::ln_one();
            let mut prob_alt = LogProb::ln_one();
            for (base, pos) in bases.into_iter().zip(start..variant.end(start)) {
                if let Some(qpos) = cigar.read_pos(pos, false, false)? {
                    let read_base = record.seq()[qpos as usize];
                    let base_qual = record.qual()[qpos as usize];
                    prob_alt += prob_read_base(read_base, *base, base_qual);
                    prob_ref += prob_read_base(read_base, ref_seq[pos as usize], base_qual);
                } else {
                    // a read that spans an SNV might have the respective position deleted (Cigar op 'D')
                    // or reference skipped (Cigar op 'N'), and the library should not choke on those reads
                    // but instead needs to know NOT to add those reads (as observations) further up
                    return Ok(None);
                }
            }
            Ok(Some((prob_ref, prob_alt)))
        } else {
            panic!("bug: unsupported variant");
        }
    }

    fn prob_sample_alt(&self, _: u32, _: &Variant, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

/// Width of band around alignment with optimal edit distance.
pub const EDIT_BAND: usize = 2;

/// Calculate read evidence for an indel.
#[derive(Debug)]
pub struct IndelEvidence {
    gap_params: IndelGapParams,
    pairhmm: pairhmm::PairHMM,
    max_window: u32,
}

impl IndelEvidence {
    /// Create a new instance.
    pub fn new(
        prob_insertion_artifact: LogProb,
        prob_deletion_artifact: LogProb,
        prob_insertion_extend_artifact: LogProb,
        prob_deletion_extend_artifact: LogProb,
        max_window: u32,
    ) -> Self {
        IndelEvidence {
            gap_params: IndelGapParams {
                prob_insertion_artifact: prob_insertion_artifact,
                prob_deletion_artifact: prob_deletion_artifact,
                prob_insertion_extend_artifact: prob_insertion_extend_artifact,
                prob_deletion_extend_artifact: prob_deletion_extend_artifact,
            },
            pairhmm: pairhmm::PairHMM::new(),
            max_window,
        }
    }

    /// Calculate probability of a certain allele.
    fn prob_allele<E: pairhmm::EmissionParameters + RefBaseEmission>(
        &mut self,
        mut allele_params: E,
        edit_dist: &EditDistanceCalculation,
    ) -> LogProb {
        let hit = edit_dist.calc_best_hit(&allele_params);
        if hit.dist == 0 {
            // METHOD: In case of a perfect match, we just take the base quality product.
            // All alternative paths in the HMM will anyway be much worse.
            allele_params.read_emission().certainty_est()
        } else {
            // METHOD: We shrink the area to run the HMM against to an environment around the best
            // edit distance hits.
            allele_params.shrink_to_hit(&hit);

            // METHOD: Further, we run the HMM on a band around the best edit distance.
            self.pairhmm.prob_related(
                &self.gap_params,
                &allele_params,
                Some(hit.dist_upper_bound()),
            )
        }
    }
}

impl AbstractReadEvidence for IndelEvidence {
    /// Calculate probability for reference and alternative indel allele. Always returns Some().
    fn prob(
        &mut self,
        record: &bam::Record,
        cigar: &CigarStringView,
        start: u32,
        variant: &Variant,
        ref_seq: &[u8],
    ) -> Result<Option<(LogProb, LogProb)>, Box<dyn Error>> {
        let read_seq = record.seq();
        let read_qual = record.qual();

        let (read_offset, read_end, breakpoint, overlap) = {
            let (varstart, varend) = match variant {
                &Variant::Deletion(_) => (start, start + variant.len()),
                &Variant::Insertion(_) => (start, start + 1),
                //TODO: add support for &Variant::Ref if we want to check against potential indel alt alleles
                &Variant::SNV(_) | &Variant::MNV(_) | &Variant::None => panic!("bug: unsupported variant"),
            };

            match (
                cigar.read_pos(varstart, true, true)?,
                cigar.read_pos(varend, true, true)?,
            ) {
                // read encloses variant
                (Some(qstart), Some(qend)) => {
                    let qstart = qstart as usize;
                    // exclusive end of variant
                    let qend = qend as usize;
                    // ensure that distance between qstart and qend does not make the window too
                    // large
                    let max_window = (self.max_window as usize).saturating_sub((qend - qstart) / 2);
                    let mut read_offset = qstart.saturating_sub(max_window);
                    let mut read_end = cmp::min(qend + max_window as usize, read_seq.len());

                    // correct for reads that enclose the entire variant while that exceeds the maximum pattern len
                    let exceed = (read_end - read_offset)
                        .saturating_sub(EditDistanceCalculation::max_pattern_len());
                    if exceed > 0 {
                        read_offset += exceed / 2;
                        read_end -= (exceed as f64 / 2.0).ceil() as usize;
                    }
                    (read_offset, read_end, varstart as usize, true)
                }

                (Some(qstart), None) => {
                    let qstart = qstart as usize;
                    let read_offset = qstart.saturating_sub(self.max_window as usize);
                    let read_end = cmp::min(qstart + self.max_window as usize, read_seq.len());
                    (read_offset, read_end, varstart as usize, true)
                }
                (None, Some(qend)) => {
                    let qend = qend as usize;
                    let read_offset = qend.saturating_sub(self.max_window as usize);
                    let read_end = cmp::min(qend + self.max_window as usize, read_seq.len());
                    (read_offset, read_end, varend as usize, true)
                }
                (None, None) => {
                    let m = read_seq.len() / 2;
                    let read_offset = m.saturating_sub(self.max_window as usize);
                    let read_end = cmp::min(m + self.max_window as usize - 1, read_seq.len());
                    let breakpoint = record.pos() as usize + m;
                    // The following should only happen with deletions.
                    // It occurs if the read comes from ref allele and is mapped within start
                    // and end of deletion. Usually, such reads strongly support the ref allele.
                    let read_enclosed_by_variant =
                        record.pos() >= varstart as i32 && cigar.end_pos() <= varend as i32;
                    (read_offset, read_end, breakpoint, read_enclosed_by_variant)
                }
            }
        };

        let start = start as usize;
        // the window on the reference should be a bit larger to allow some flexibility with close
        // indels. But it should not be so large that the read can align outside of the breakpoint.
        let ref_window = (self.max_window as f64 * 1.5) as usize;

        // read emission
        let read_emission = ReadEmission::new(&read_seq, read_qual, read_offset, read_end);
        let edit_dist = EditDistanceCalculation::new((read_offset..read_end).map(|i| read_seq[i]));

        if !overlap {
            // If there is no overlap, normalization below would anyway lead to 0.5 vs 0.5,
            // multiplied with certainty estimate. Hence, we can skip the entire HMM calculation!
            let p = LogProb::from(Prob(0.5));
            return Ok(Some((p, p)));
        }

        // ref allele
        let mut prob_ref = self.prob_allele(
            ReferenceEmissionParams {
                ref_seq: ref_seq,
                ref_offset: breakpoint.saturating_sub(ref_window),
                ref_end: cmp::min(breakpoint + ref_window, ref_seq.len()),
                read_emission: &read_emission,
            },
            &edit_dist,
        );

        // alt allele
        let mut prob_alt = match variant {
            &Variant::Deletion(_) => self.prob_allele(
                DeletionEmissionParams {
                    ref_seq: ref_seq,
                    ref_offset: start.saturating_sub(ref_window),
                    ref_end: cmp::min(start + ref_window, ref_seq.len() - variant.len() as usize),
                    del_start: start,
                    del_len: variant.len() as usize,
                    read_emission: &read_emission,
                },
                &edit_dist,
            ),
            &Variant::Insertion(ref ins_seq) => {
                let l = ins_seq.len() as usize;

                self.prob_allele(
                    InsertionEmissionParams {
                        ref_seq: ref_seq,
                        ref_offset: start.saturating_sub(ref_window),
                        ref_end: cmp::min(start + l + ref_window, ref_seq.len()),
                        ins_start: start,
                        ins_len: l,
                        ins_end: start + l,
                        ins_seq: ins_seq,
                        read_emission: &read_emission,
                    },
                    &edit_dist,
                )
            }
            &Variant::SNV(_) | &Variant::MNV(_) | &Variant::None => {
                panic!("bug: unsupported variant");
            }
        };
        assert!(!prob_ref.is_nan());
        assert!(!prob_alt.is_nan());

        // METHOD: Normalize probabilities. By this, we avoid biases due to proximal variants that are in
        // cis with the considered one. They are normalized away since they affect both ref and alt.
        // In a sense, this assumes that the two considered alleles are the only possible ones.
        // However, if the read actually comes from a third allele, both probabilities will be
        // equally bad, and the normalized one will not prefer any of them.
        // This is ok, because for the likelihood function only the ratio between the two
        // probabilities is relevant!

        if prob_ref != LogProb::ln_zero() && prob_alt != LogProb::ln_zero() {
            // METHOD: Only perform normalization if both probs are non-zero
            // otherwise, we would artificially magnify the ratio
            // (compared to an epsilon for the zero case).
            let prob_total = prob_alt.ln_add_exp(prob_ref);
            prob_ref -= prob_total;
            prob_alt -= prob_total;
        }

        if prob_ref == LogProb::ln_zero() && prob_alt == LogProb::ln_zero() {
            // METHOD: if both are zero, use 0.5 instead. Since only the ratio matters, this
            // has the same effect, without rendering the entire pileup likelihood zero.
            prob_ref = LogProb::from(Prob(0.5));
            prob_alt = prob_ref;
        }

        Ok(Some((prob_ref, prob_alt)))
    }

    /// Probability to sample read from alt allele given the average feasible positions observed
    /// from a subsample of the mapped reads.
    ///
    /// The key idea is calculate the probability as number of valid placements (considering the
    /// max softclip allowed by the mapper) over all possible placements.
    fn prob_sample_alt(
        &self,
        read_len: u32,
        variant: &Variant,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        // TODO for long reads, always return One
        let delta = match variant {
            &Variant::Deletion(_) => variant.len() as u32,
            &Variant::Insertion(_) => variant.len() as u32,
            &Variant::SNV(_) | &Variant::MNV(_) | &Variant::None => panic!("unsupported variant"),
        };

        let feasible = alignment_properties.feasible_bases(read_len, variant);

        let prob = {
            let n_alt = cmp::min(delta, read_len);
            let n_alt_valid = cmp::min(n_alt, feasible);

            LogProb((n_alt_valid as f64).ln() - (n_alt as f64).ln())
        };
        assert!(prob.is_valid());

        prob
    }
}

lazy_static! {
    static ref PROB_CONFUSION: LogProb = LogProb::from(Prob(0.3333));
}

/// Calculate probability of read_base given ref_base.
pub fn prob_read_base(read_base: u8, ref_base: u8, base_qual: u8) -> LogProb {
    let prob_miscall = prob_read_base_miscall(base_qual);

    if read_base.to_ascii_uppercase() == ref_base.to_ascii_uppercase() {
        prob_miscall.ln_one_minus_exp()
    } else {
        // TODO replace the second term with technology specific confusion matrix
        prob_miscall + *PROB_CONFUSION
    }
}

/// unpack miscall probability of read_base.
pub fn prob_read_base_miscall(base_qual: u8) -> LogProb {
    LogProb::from(PHREDProb::from((base_qual) as f64))
}

/// Gap parameters for PairHMM.
#[derive(Debug, Clone)]
pub struct IndelGapParams {
    pub prob_insertion_artifact: LogProb,
    pub prob_deletion_artifact: LogProb,
    pub prob_insertion_extend_artifact: LogProb,
    pub prob_deletion_extend_artifact: LogProb,
}

impl pairhmm::GapParameters for IndelGapParams {
    #[inline]
    fn prob_gap_x(&self) -> LogProb {
        self.prob_insertion_artifact
    }

    #[inline]
    fn prob_gap_y(&self) -> LogProb {
        self.prob_deletion_artifact
    }

    #[inline]
    fn prob_gap_x_extend(&self) -> LogProb {
        self.prob_insertion_extend_artifact
    }

    #[inline]
    fn prob_gap_y_extend(&self) -> LogProb {
        self.prob_deletion_extend_artifact
    }
}

impl pairhmm::StartEndGapParameters for IndelGapParams {
    /// Semiglobal alignment: return true.
    #[inline]
    fn free_start_gap_x(&self) -> bool {
        true
    }

    /// Semiglobal alignment: return true.
    #[inline]
    fn free_end_gap_x(&self) -> bool {
        true
    }

    /// Semiglobal alignment: return 1.0.
    #[inline]
    fn prob_start_gap_x(&self, _: usize) -> LogProb {
        LogProb::ln_one()
    }
}

macro_rules! default_emission {
    () => (
        #[inline]
        fn prob_emit_xy(&self, i: usize, j: usize) -> pairhmm::XYEmission {
            let r = self.ref_base(i);
            self.read_emission.prob_match_mismatch(j, r)
        }

        #[inline]
        fn prob_emit_x(&self, _: usize) -> LogProb {
            LogProb::ln_one()
        }

        #[inline]
        fn prob_emit_y(&self, j: usize) -> LogProb {
            self.read_emission.prob_insertion(j)
        }

        #[inline]
        fn len_y(&self) -> usize {
            self.read_emission.read_end - self.read_emission.read_offset
        }
    )
}

#[derive(Debug)]
pub struct ReadEmission<'a> {
    read_seq: &'a bam::record::Seq<'a>,
    any_miscall: Vec<LogProb>,
    no_miscall: Vec<LogProb>,
    read_offset: usize,
    read_end: usize,
}

impl<'a> ReadEmission<'a> {
    pub fn new(
        read_seq: &'a bam::record::Seq<'a>,
        qual: &[u8],
        read_offset: usize,
        read_end: usize,
    ) -> Self {
        let mut any_miscall = vec![LogProb::ln_zero(); read_end - read_offset];
        let mut no_miscall = any_miscall.clone();
        for (j, j_) in (read_offset..read_end).enumerate() {
            let prob_miscall = prob_read_base_miscall(qual[j_]);
            any_miscall[j] = prob_miscall;
            no_miscall[j] = prob_miscall.ln_one_minus_exp();
        }
        ReadEmission {
            read_seq,
            any_miscall,
            no_miscall,
            read_offset,
            read_end,
        }
    }

    fn particular_miscall(&self, j: usize) -> LogProb {
        self.any_miscall[j] + *PROB_CONFUSION
    }

    /// Calculate probability of read_base given ref_base.
    pub fn prob_match_mismatch(&self, j: usize, ref_base: u8) -> pairhmm::XYEmission {
        let read_base = self.read_seq[self.project_j(j)];

        if read_base.to_ascii_uppercase() == ref_base.to_ascii_uppercase() {
            pairhmm::XYEmission::Match(self.no_miscall[j])
        } else {
            // TODO replace the second term with technology specific confusion matrix
            pairhmm::XYEmission::Mismatch(self.particular_miscall(j))
        }
    }

    pub fn prob_insertion(&self, j: usize) -> LogProb {
        self.any_miscall[j]
    }

    #[inline]
    fn project_j(&self, j: usize) -> usize {
        j + self.read_offset
    }

    /// Calculate probability that none of the bases is miscalled.
    pub fn certainty_est(&self) -> LogProb {
        self.no_miscall.iter().sum()
    }
}

pub trait RefBaseEmission: Debug {
    fn ref_base(&self, i: usize) -> u8;

    fn ref_offset(&self) -> usize;

    fn ref_end(&self) -> usize;

    fn set_ref_offset(&mut self, value: usize);

    fn set_ref_end(&mut self, value: usize);

    fn read_emission(&self) -> &ReadEmission;

    fn shrink_to_hit(&mut self, hit: &EditDistanceHit) {
        self.set_ref_end(cmp::min(
            self.ref_offset() + hit.end + EDIT_BAND,
            self.ref_end(),
        ));
        self.set_ref_offset(self.ref_offset() + hit.start.saturating_sub(EDIT_BAND));
    }
}

macro_rules! default_ref_base_emission {
    () => (
        fn ref_offset(&self) -> usize {
            self.ref_offset
        }

        fn ref_end(&self) -> usize {
            self.ref_end
        }

        fn set_ref_offset(&mut self, value: usize) {
            self.ref_offset = value;
        }

        fn set_ref_end(&mut self, value: usize) {
            self.ref_end = value;
        }

        fn read_emission(&self) -> &ReadEmission {
            self.read_emission
        }
    )
}

/// Emission parameters for PairHMM over reference allele.
#[derive(Debug)]
pub struct ReferenceEmissionParams<'a> {
    ref_seq: &'a [u8],
    ref_offset: usize,
    ref_end: usize,
    read_emission: &'a ReadEmission<'a>,
}

impl<'a> RefBaseEmission for ReferenceEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        self.ref_seq[i + self.ref_offset]
    }

    default_ref_base_emission!();
}

impl<'a> pairhmm::EmissionParameters for ReferenceEmissionParams<'a> {
    default_emission!();

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
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

/// Emission parameters for PairHMM over insertion allele.
#[derive(Debug)]
pub struct InsertionEmissionParams<'a> {
    ref_seq: &'a [u8],
    ref_offset: usize,
    ref_end: usize,
    ins_start: usize,
    ins_end: usize,
    ins_len: usize,
    ins_seq: &'a [u8],
    read_emission: &'a ReadEmission<'a>,
}

impl<'a> RefBaseEmission for InsertionEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;
        if i_ <= self.ins_start {
            self.ref_seq[i_]
        } else if i_ > self.ins_end {
            self.ref_seq[i_ - self.ins_len]
        } else {
            self.ins_seq[i_ - (self.ins_start + 1)]
        }
    }

    default_ref_base_emission!();
}

impl<'a> pairhmm::EmissionParameters for InsertionEmissionParams<'a> {
    default_emission!();

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset + self.ins_len
    }
}

pub struct EditDistanceCalculation {
    myers: Myers<u128>,
    read_seq_len: usize,
}

impl EditDistanceCalculation {
    pub fn max_pattern_len() -> usize {
        128
    }

    /// Create new instance.
    ///
    /// # Arguments
    /// * `read_seq` - read sequence in window (may not exceed 128 bases).
    pub fn new<P>(read_seq: P) -> Self
    where
        P: Iterator<Item = u8> + DoubleEndedIterator + ExactSizeIterator,
    {
        let l = read_seq.len();
        EditDistanceCalculation {
            myers: Myers::new(read_seq.rev()),
            read_seq_len: l,
        }
    }

    /// Returns a reasonable upper bound for the edit distance in order to band the pairHMM computation.
    /// We use the best edit distance and add 5.
    pub fn calc_best_hit<E: pairhmm::EmissionParameters + RefBaseEmission>(
        &self,
        emission_params: &E,
    ) -> EditDistanceHit {
        let ref_seq = (0..emission_params.len_x())
            .rev()
            .map(|i| emission_params.ref_base(i).to_ascii_uppercase());
        let mut best_dist = u8::max_value();
        let mut positions = Vec::new();
        for (pos, dist) in self.myers.find_all_end(ref_seq, u8::max_value()) {
            if dist < best_dist {
                positions.clear();
                positions.push(pos);
                best_dist = dist;
            } else if dist == best_dist {
                positions.push(pos);
            }
        }
        let ambiguous = positions.len() > 1;

        // We find a pos relative to ref end, hence we have to project it to a position relative to
        // the start.
        let project = |pos| emission_params.len_x() - pos;
        let start = project(*positions.last().unwrap()).saturating_sub(best_dist as usize);
        // take the last (aka first because we are mapping backwards) position for an upper bound of the putative end
        let end = cmp::min(
            project(positions[0]) + self.read_seq_len + best_dist as usize,
            emission_params.len_x(),
        );
        EditDistanceHit {
            start,
            end,
            ambiguous,
            dist: best_dist,
        }
    }
}

#[derive(Debug, Clone)]
pub struct EditDistanceHit {
    start: usize,
    end: usize,
    dist: u8,
    ambiguous: bool,
}

impl EditDistanceHit {
    pub fn dist_upper_bound(&self) -> usize {
        self.dist as usize + EDIT_BAND
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::model;

    use rust_htslib::bam::record::{Cigar, CigarString};
    use std::str;

    #[test]
    fn test_prob_none() {
        let ref_seq: Vec<u8> = b"GATTACA"[..].to_owned();

        let mut records: Vec<bam::Record> = Vec::new();
        let mut qname: &[u8];
        let mut seq: &[u8];

        let mut none_evidence = NoneEvidence::new();

        // Ignore leading HardClip, skip leading SoftClip, reference nucleotide
        qname = b"HC_SC_ref";
        let cigar = CigarString(vec![
            Cigar::HardClip(3),
            Cigar::SoftClip(1),
            Cigar::Match(5),
        ]);
        seq = b"TATTaC";
        let qual = [20, 30, 30, 30, 40, 30];
        let mut record1 = bam::Record::new();
        record1.set(qname, Some(&cigar), seq, &qual);
        record1.set_pos(1);
        records.push(record1);

        // Ignore leading HardClip, skip leading SoftClip, non-reference nucleotide
        qname = b"HC_SC_non-ref";
        let cigar = CigarString(vec![
            Cigar::HardClip(5),
            Cigar::SoftClip(2),
            Cigar::Match(4),
        ]);
        seq = b"TTTTCC";
        let qual = [15, 15, 20, 20, 30, 20];
        let mut record2 = bam::Record::new();
        record2.set(qname, Some(&cigar), seq, &qual);
        record2.set_pos(2);
        records.push(record2);

        // reference nucleotide, trailing SoftClip, trailing HardClip
        qname = b"ref_SC_HC";
        let cigar = CigarString(vec![
            Cigar::Match(3),
            Cigar::SoftClip(2),
            Cigar::HardClip(7),
        ]);
        seq = b"ACATA";
        let qual = [50, 20, 20, 20, 20, 20];
        let mut record3 = bam::Record::new();
        record3.set(qname, Some(&cigar), seq, &qual);
        record3.set_pos(4);
        records.push(record3);

        // three nucleotide Deletion covering Ref position
        qname = b"M_3Del_M";
        let cigar = CigarString(vec![Cigar::Match(3), Cigar::Del(3), Cigar::Match(1)]);
        seq = b"GATA";
        let qual = [10, 30, 30, 30];
        let mut record4 = bam::Record::new();
        record4.set(qname, Some(&cigar), seq, &qual);
        record4.set_pos(0);
        records.push(record4);

        // truth
        let probs_ref = [0.9999, 0.001, 0.99999];
        let probs_alt = [0.0001, 0.999, 0.00001];
        let eps = [0.00001, 0.0001, 0.000001];

        let vpos = 4;
        let variant = model::Variant::None;
        for (i, mut rec) in records.into_iter().enumerate() {
            rec.cache_cigar();
            println!("{}", str::from_utf8(rec.qname()).unwrap());
            if let Ok(Some((prob_ref, prob_alt))) =
                none_evidence.prob(&rec, rec.cigar_cached().unwrap(), vpos, &variant, &ref_seq)
            {
                println!("{:?}", rec.cigar_cached());
                println!(
                    "Pr(ref)={} Pr(alt)={}",
                    (*prob_ref).exp(),
                    (*prob_alt).exp()
                );
                assert_relative_eq!((*prob_ref).exp(), probs_ref[i], epsilon = eps[i]);
                assert_relative_eq!((*prob_alt).exp(), probs_alt[i], epsilon = eps[i]);
            } else {
                // tests for reference position not being covered should be pushed onto records last
                // and should have 10 as the quality value of the first base in seq
                assert_eq!(rec.qual()[0], 10);
            }
        }
    }
    #[test]
    fn test_prob_snv() {
        let ref_seq: Vec<u8> = b"CCTATACGCGT"[..].to_owned();

        let mut records: Vec<bam::Record> = Vec::new();
        let mut qname: &[u8];
        let mut seq: &[u8];

        let mut snv_evidence = SNVEvidence::new();

        // Ignore leading HardClip, skip leading SoftClip, reference nucleotide
        qname = b"HC_SC_M";
        let cigar = CigarString(vec![
            Cigar::HardClip(5),
            Cigar::SoftClip(2),
            Cigar::Match(6),
        ]);
        seq = b"AATATACG";
        let qual = [20, 20, 30, 30, 30, 40, 30, 30];
        let mut record1 = bam::Record::new();
        record1.set(qname, Some(&cigar), seq, &qual);
        record1.set_pos(2);
        records.push(record1);

        // Ignore leading HardClip, skip leading Insertion, alternative nucleotide
        qname = b"HC_Ins_M";
        let cigar = CigarString(vec![Cigar::HardClip(2), Cigar::Ins(2), Cigar::Match(6)]);
        seq = b"TTTATGCG";
        let qual = [20, 20, 20, 20, 20, 30, 20, 20];
        let mut record2 = bam::Record::new();
        record2.set(qname, Some(&cigar), seq, &qual);
        record2.set_pos(2);
        records.push(record2);

        // Matches and deletion before position, reference nucleotide
        qname = b"Eq_Diff_Del_Eq";
        let cigar = CigarString(vec![
            Cigar::Equal(2),
            Cigar::Diff(1),
            Cigar::Del(2),
            Cigar::Equal(5),
        ]);
        seq = b"CCAACGCG";
        let qual = [30, 30, 30, 50, 30, 30, 30, 30];
        let mut record3 = bam::Record::new();
        record3.set(qname, Some(&cigar), seq, &qual);
        record3.set_pos(0);
        records.push(record3);

        // single nucleotide Deletion covering SNV position
        qname = b"M_Del_M";
        let cigar = CigarString(vec![Cigar::Match(4), Cigar::Del(1), Cigar::Match(4)]);
        seq = b"CTATCGCG";
        let qual = [10, 30, 30, 30, 30, 30, 30, 30];
        let mut record4 = bam::Record::new();
        record4.set(qname, Some(&cigar), seq, &qual);
        record4.set_pos(1);
        records.push(record4);

        // three nucleotide RefSkip covering SNV position
        qname = b"M_RefSkip_M";
        let cigar = CigarString(vec![
            Cigar::Equal(1),
            Cigar::Diff(1),
            Cigar::Equal(2),
            Cigar::RefSkip(3),
            Cigar::Match(4),
        ]);
        seq = b"CTTAGCGT";
        let qual = [10, 30, 30, 30, 30, 30, 30, 30];
        let mut record5 = bam::Record::new();
        record5.set(qname, Some(&cigar), seq, &qual);
        record5.set_pos(0);
        records.push(record5);

        // truth
        let probs_ref = [0.9999, 0.00033, 0.99999];
        let probs_alt = [0.000033, 0.999, 0.0000033];
        let eps = [0.000001, 0.00001, 0.0000001];

        let vpos = 5;
        let variant = model::Variant::SNV(b'G');
        for (i, mut rec) in records.into_iter().enumerate() {
            rec.cache_cigar();
            println!("{}", str::from_utf8(rec.qname()).unwrap());
            if let Ok(Some((prob_ref, prob_alt))) =
                snv_evidence.prob(&rec, rec.cigar_cached().unwrap(), vpos, &variant, &ref_seq)
            {
                println!("{:?}", rec.cigar_cached());
                println!(
                    "Pr(ref)={} Pr(alt)={}",
                    (*prob_ref).exp(),
                    (*prob_alt).exp()
                );
                assert_relative_eq!((*prob_ref).exp(), probs_ref[i], epsilon = eps[i]);
                assert_relative_eq!((*prob_alt).exp(), probs_alt[i], epsilon = eps[i]);
            } else {
                // tests for reference position not being covered should be pushed onto records last
                // and should have 10 as the quality value of the first base in seq
                assert_eq!(rec.qual()[0], 10);
            }
        }
    }
}

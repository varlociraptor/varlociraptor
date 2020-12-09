// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::cmp::Ordering;
use std::ops::Range;
use std::rc::Rc;
use std::str;
use std::sync::Arc;
use std::usize;

use anyhow::Result;
use bio::alignment::AlignmentOperation;
use bio::stats::{self, pairhmm::GapParameters, pairhmm::PairHMM, LogProb, Prob};
use bio_types::genome;
use bio_types::genome::AbstractInterval;
use rust_htslib::bam;

use crate::errors::Error;
use crate::reference;
use crate::utils;
use crate::variants::evidence::observation::Strand;
use crate::variants::evidence::realignment::edit_distance::EditDistanceCalculation;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, ReferenceEmissionParams};
use crate::variants::types::{AlleleSupport, AlleleSupportBuilder, SingleLocus};

pub(crate) mod edit_distance;
pub(crate) mod pairhmm;

use crate::variants::evidence::realignment::edit_distance::EditDistanceHit;

pub(crate) struct CandidateRegion {
    overlap: bool,
    read_interval: Range<usize>,
    ref_interval: Range<usize>,
}

pub(crate) trait Realignable<'a> {
    type EmissionParams: stats::pairhmm::EmissionParameters + pairhmm::RefBaseEmission;

    fn alt_emission_params(
        &self,
        read_emission_params: Rc<ReadEmission<'a>>,
        ref_buffer: Arc<reference::Buffer>,
        ref_interval: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<Self::EmissionParams>>;

    /// Returns true if reads emitted from alt allele
    /// may be interpreted as revcomp reads by the mapper.
    /// In such a case, the realigner needs to consider both
    /// forward and reverse sequence of the read.
    fn maybe_revcomp(&self) -> bool {
        false
    }
}

pub(crate) trait Realigner {
    fn candidate_region(
        &self,
        record: &bam::Record,
        locus: &genome::Interval,
    ) -> Result<CandidateRegion> {
        let cigar = record.cigar_cached().unwrap();

        let locus_start = locus.range().start;
        let locus_end = locus.range().end;

        let ref_seq = self.ref_buffer().seq(locus.contig())?;

        let ref_interval = |breakpoint: usize| {
            breakpoint.saturating_sub(self.ref_window())
                ..cmp::min(breakpoint + self.ref_window(), ref_seq.len())
        };

        Ok(
            match (
                cigar.read_pos(locus_start as u32, true, true)?,
                cigar.read_pos(locus_end as u32, true, true)?,
            ) {
                // read encloses variant
                (Some(qstart), Some(qend)) => {
                    let qstart = qstart as usize;
                    // exclusive end of variant
                    let qend = qend as usize;
                    // ensure that distance between qstart and qend does not make the window too
                    // large
                    let max_window =
                        (self.max_window() as usize).saturating_sub((qend - qstart) / 2);
                    let mut read_offset = qstart.saturating_sub(max_window);
                    let mut read_end = cmp::min(qend + max_window as usize, record.seq_len());

                    // correct for reads that enclose the entire variant while that exceeds the maximum pattern len
                    let exceed = (read_end - read_offset)
                        .saturating_sub(EditDistanceCalculation::max_pattern_len());
                    if exceed > 0 {
                        read_offset += exceed / 2;
                        read_end -= (exceed as f64 / 2.0).ceil() as usize;
                    }

                    CandidateRegion {
                        overlap: true,
                        read_interval: read_offset..read_end,
                        ref_interval: ref_interval(locus_start as usize),
                    }
                }

                // read overlaps from right
                (Some(qstart), None) => {
                    let qstart = qstart as usize;
                    let read_offset = qstart.saturating_sub(self.max_window() as usize);
                    let read_end = cmp::min(qstart + self.max_window() as usize, record.seq_len());

                    CandidateRegion {
                        overlap: true,
                        read_interval: read_offset..read_end,
                        ref_interval: ref_interval(locus_start as usize),
                    }
                }

                // read overlaps from left
                (None, Some(qend)) => {
                    let qend = qend as usize;
                    let read_offset = qend.saturating_sub(self.max_window() as usize);
                    let read_end = cmp::min(qend + self.max_window() as usize, record.seq_len());

                    CandidateRegion {
                        overlap: true,
                        read_interval: read_offset..read_end,
                        ref_interval: ref_interval(locus_end as usize),
                    }
                }

                // no overlap
                (None, None) => {
                    let m = record.seq_len() / 2;
                    let read_offset = m.saturating_sub(self.max_window() as usize);
                    let read_end = cmp::min(m + self.max_window() as usize - 1, record.seq_len());
                    let breakpoint = record.pos() as usize + m;
                    // The following should only happen with deletions.
                    // It occurs if the read comes from ref allele and is mapped within start
                    // and end of deletion. Usually, such reads strongly support the ref allele.
                    let read_enclosed_by_variant =
                        record.pos() >= locus_start as i64 && cigar.end_pos() <= locus_end as i64;

                    CandidateRegion {
                        overlap: read_enclosed_by_variant,
                        read_interval: read_offset..read_end,
                        ref_interval: ref_interval(breakpoint),
                    }
                }
            },
        )
    }

    fn ref_window(&self) -> usize {
        // METHOD: the window on the reference should be a bit larger to allow some flexibility with close
        // indels. But it should not be so large that the read can align outside of the breakpoint.
        (self.max_window() as f64 * 1.5) as usize
    }

    fn allele_support<'a, V, L>(
        &mut self,
        record: &'a bam::Record,
        loci: L,
        variant: &V,
    ) -> Result<AlleleSupport>
    where
        V: Realignable<'a>,
        L: IntoIterator,
        L::Item: AsRef<SingleLocus>,
    {
        // Obtain candidate regions from matching loci.
        let candidate_regions: Result<Vec<_>> = loci
            .into_iter()
            .filter_map(|locus| {
                if locus.as_ref().contig() == record.contig() {
                    Some(self.candidate_region(record, locus.as_ref()))
                } else {
                    None
                }
            })
            .collect();
        let candidate_regions = candidate_regions?;

        // Check if anything overlaps.
        if !candidate_regions.iter().any(|region| region.overlap) {
            // If there is no overlap, normalization below would anyway lead to 0.5 vs 0.5,
            // multiplied with certainty estimate. Hence, we can skip the entire HMM calculation!
            let p = LogProb(0.5f64.ln());
            return Ok(AlleleSupportBuilder::default()
                .prob_ref_allele(p)
                .prob_alt_allele(p)
                .strand(Strand::None)
                .build()
                .unwrap());
        }

        // Merge overlapping reference regions of all read-overlapping candidate regions.
        // Regions are sorted by reference position, hence we can merge from left to right.
        let mut merged_regions = Vec::new();
        for region in candidate_regions
            .into_iter()
            .filter(|region| region.overlap)
        {
            if merged_regions.is_empty() {
                merged_regions.push(region);
            } else {
                let last = merged_regions.last_mut().unwrap();
                if region.ref_interval.start <= last.ref_interval.end {
                    // They overlap, hence merge.
                    last.ref_interval = last.ref_interval.start..region.ref_interval.end;
                    last.read_interval =
                        cmp::min(last.read_interval.start, region.read_interval.start)
                            ..cmp::max(last.read_interval.end, region.read_interval.end)
                } else {
                    // No overlap, hence push.
                    merged_regions.push(region);
                }
            }
        }

        // Calculate independent probabilities over all merged regions.
        let mut prob_ref_all = LogProb::ln_one();
        let mut prob_alt_all = LogProb::ln_one();

        let ref_seq = self.ref_buffer().seq(record.contig())?;
        let read_seq: bam::record::Seq<'a> = record.seq();
        let read_qual = record.qual();

        let aux_strand_info = utils::aux_tag_strand_info(record);
        let mut strand = Strand::None;

        for region in merged_regions {
            // read emission
            let read_emission = Rc::new(ReadEmission::new(
                read_seq,
                read_qual,
                region.read_interval.start,
                region.read_interval.end,
            ));
            let mut edit_dist =
                EditDistanceCalculation::new(region.read_interval.clone().map(|i| read_seq[i]));

            // ref allele
            let mut prob_ref = self.prob_allele(
                &mut [ReferenceEmissionParams {
                    ref_seq: Arc::clone(&ref_seq),
                    ref_offset: region.ref_interval.start,
                    ref_end: region.ref_interval.end,
                    read_emission: Rc::clone(&read_emission),
                }],
                &mut edit_dist,
            );

            let mut prob_alt = self.prob_allele(
                &mut variant.alt_emission_params(
                    Rc::clone(&read_emission),
                    Arc::clone(self.ref_buffer()),
                    &genome::Interval::new(
                        record.contig().to_owned(),
                        region.ref_interval.start as u64..region.ref_interval.end as u64,
                    ),
                    self.ref_window(),
                )?,
                &mut edit_dist,
            );

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
                debug!(
                    "Record {} neither supports reference nor alternative allele during realignment.",
                    str::from_utf8(record.qname()).unwrap()
                );
            }

            if prob_ref != prob_alt {
                // METHOD: if record is not informative, we don't want to
                // retain its information (e.g. strand).
                if let Some(strand_info) = aux_strand_info {
                    if let Some(s) = strand_info.get(region.read_interval) {
                        strand |= Strand::from_aux(s)?;
                    } else {
                        return Err(Error::ReadPosOutOfBounds.into());
                    }
                }
            }

            // METHOD: probabilities of independent regions are combined here.
            prob_ref_all += prob_ref;
            prob_alt_all += prob_alt;
        }

        if aux_strand_info.is_none() && prob_ref_all != prob_alt_all {
            // METHOD: if record is not informative, we don't want to
            // retain its information (e.g. strand).
            strand = Strand::from_record(record);
        }

        Ok(AlleleSupportBuilder::default()
            .strand(strand)
            .prob_ref_allele(prob_ref_all)
            .prob_alt_allele(prob_alt_all)
            .build()
            .unwrap())
    }

    /// Calculate probability of a certain allele.
    fn prob_allele<E>(
        &mut self,
        candidate_allele_params: &mut [E],
        edit_dist: &mut edit_distance::EditDistanceCalculation,
    ) -> LogProb
    where
        E: stats::pairhmm::EmissionParameters + pairhmm::RefBaseEmission,
    {
        let mut hits = Vec::new();
        let mut best_dist = None;
        for params in candidate_allele_params.iter_mut() {
            if let Some(hit) = edit_dist.calc_best_hit(params, None) {
                match best_dist.map_or(Ordering::Less, |best_dist| hit.dist().cmp(&best_dist)) {
                    Ordering::Less => {
                        hits.clear();
                        best_dist = Some(hit.dist());
                        hits.push((hit, params));
                    }
                    Ordering::Equal => {
                        hits.push((hit, params));
                    }
                    Ordering::Greater => (),
                }
            }
        }

        let mut prob = None;
        // METHOD: for equal best edit dists, we have to compare the probabilities and take the best.
        for (hit, allele_params) in hits {
            if hit.dist() == 0 {
                // METHOD: In case of a perfect match, we just take the base quality product.
                // All alternative paths in the HMM will anyway be much worse.
                prob = Some(allele_params.read_emission().certainty_est());
            } else {
                let p = self.calculate_prob_allele(&hit, allele_params);

                if prob.map_or(true, |prob| p > prob) {
                    prob.replace(p);
                }
            }
        }

        // This is safe, as there will be always one probability at least.
        prob.unwrap()
    }

    fn calculate_prob_allele<E>(&mut self, hit: &EditDistanceHit, allele_params: &mut E) -> LogProb
    where
        E: stats::pairhmm::EmissionParameters + pairhmm::RefBaseEmission;

    fn ref_buffer(&self) -> &Arc<reference::Buffer>;

    fn max_window(&self) -> u64;
}

#[derive(Clone)]
pub(crate) struct PairHMMRealigner {
    gap_params: pairhmm::GapParams,
    pairhmm: PairHMM,
    max_window: u64,
    ref_buffer: Arc<reference::Buffer>,
}

impl PairHMMRealigner {
    /// Create a new instance.
    pub(crate) fn new(
        ref_buffer: Arc<reference::Buffer>,
        gap_params: pairhmm::GapParams,
        max_window: u64,
    ) -> Self {
        let pairhmm = PairHMM::new(&gap_params);
        PairHMMRealigner {
            gap_params,
            pairhmm,
            max_window,
            ref_buffer,
        }
    }
}

impl Realigner for PairHMMRealigner {
    fn ref_buffer(&self) -> &Arc<reference::Buffer> {
        &self.ref_buffer
    }

    fn max_window(&self) -> u64 {
        self.max_window
    }

    fn calculate_prob_allele<E>(&mut self, hit: &EditDistanceHit, allele_params: &mut E) -> LogProb
    where
        E: stats::pairhmm::EmissionParameters + pairhmm::RefBaseEmission,
    {
        // METHOD: We shrink the area to run the HMM against to an environment around the best
        // edit distance hits.
        allele_params.shrink_to_hit(&hit);

        // METHOD: Further, we run the HMM on a band around the best edit distance.
        self.pairhmm.prob_related(
            allele_params,
            &self.gap_params,
            Some(hit.dist_upper_bound()),
        )
    }
}

#[derive(Clone)]
pub(crate) struct PathHMMRealigner {
    gap_params: pairhmm::GapParams,
    max_window: u64,
    ref_buffer: Arc<reference::Buffer>,
    prob_no_gap: LogProb,
    prob_close_gap_x: LogProb,
    prob_close_gap_y: LogProb,
    prob_extend_or_reopen_gap_x: LogProb,
    prob_extend_or_reopen_gap_y: LogProb,
}

impl PathHMMRealigner {
    pub(crate) fn new(
        gap_params: pairhmm::GapParams,
        max_window: u64,
        ref_buffer: Arc<reference::Buffer>,
    ) -> Self {
        let prob_no_gap = gap_params
            .prob_gap_x()
            .ln_add_exp(gap_params.prob_gap_y())
            .ln_one_minus_exp();
        let prob_close_gap_x = gap_params.prob_gap_x_extend().ln_one_minus_exp();
        let prob_close_gap_y = gap_params.prob_gap_y_extend().ln_one_minus_exp();
        let prob_extend_or_reopen_gap_x = gap_params
            .prob_gap_x_extend()
            .ln_add_exp(prob_close_gap_x + gap_params.prob_gap_x());
        let prob_extend_or_reopen_gap_y = gap_params
            .prob_gap_y_extend()
            .ln_add_exp(prob_close_gap_y + gap_params.prob_gap_y());
        PathHMMRealigner {
            gap_params,
            max_window,
            ref_buffer,
            prob_no_gap,
            prob_close_gap_x,
            prob_close_gap_y,
            prob_extend_or_reopen_gap_x,
            prob_extend_or_reopen_gap_y,
        }
    }
}

impl Realigner for PathHMMRealigner {
    fn ref_buffer(&self) -> &Arc<reference::Buffer> {
        &self.ref_buffer
    }

    fn max_window(&self) -> u64 {
        self.max_window
    }

    fn calculate_prob_allele<E>(&mut self, hit: &EditDistanceHit, allele_params: &mut E) -> LogProb
    where
        E: stats::pairhmm::EmissionParameters + pairhmm::RefBaseEmission,
    {
        let mut best_prob = None;
        for alignment in hit.alignments() {
            let mut prob = LogProb::ln_one();
            let mut pos_ref = alignment.start();
            let mut pos_read = 0;
            let mut prev_operation = None;
            for operation in alignment.operations() {
                match operation {
                    AlignmentOperation::Match | AlignmentOperation::Subst => {
                        // transitions
                        match prev_operation {
                            Some(&AlignmentOperation::Del) => {
                                prob += self.prob_close_gap_y;
                            },
                            Some(&AlignmentOperation::Ins) => {
                                prob += self.prob_close_gap_x;
                            },
                            Some(&AlignmentOperation::Match) | Some(&AlignmentOperation::Subst) => {
                                prob += self.prob_no_gap;
                            }
                            _ => ()
                        }

                        let emission = allele_params.prob_emit_xy(pos_ref, pos_read);
                        prob += emission.prob();
                        pos_ref += 1;
                        pos_read += 1;
                    },
                    AlignmentOperation::Del => {
                        // transitions
                        match prev_operation {
                            Some(&AlignmentOperation::Del) => {
                                prob += self.prob_extend_or_reopen_gap_y;
                            },
                            Some(&AlignmentOperation::Ins) => {
                                prob += self.prob_close_gap_x + self.gap_params.prob_gap_y();
                            },
                            None | Some(&AlignmentOperation::Match) | Some(&AlignmentOperation::Subst) => {
                                prob += self.gap_params.prob_gap_y();
                            },
                            _ => (),
                        }
                        prob += allele_params.prob_emit_x(pos_ref);
                        pos_ref += 1;
                    },
                    AlignmentOperation::Ins => {
                        // transitions
                        match prev_operation {
                            Some(&AlignmentOperation::Ins) => {
                                prob += self.prob_extend_or_reopen_gap_x;
                            },
                            Some(&AlignmentOperation::Del) => {
                                prob += self.prob_close_gap_y + self.gap_params.prob_gap_x();
                            },
                            None | Some(&AlignmentOperation::Match) | Some(&AlignmentOperation::Subst) => {
                                prob += self.gap_params.prob_gap_x();
                            },
                            _ => (),
                        }
                        prob += allele_params.prob_emit_y(pos_read);
                        pos_read += 1;
                    },
                    _ => panic!("bug: unsupported alignment operation (Myers algorithm should create a semiglobal alignment)")
                }
                prev_operation = Some(operation);
            }
            // If better than best probability, keep this.
            if best_prob.map_or(true, |best| prob > best) {
                best_prob = Some(prob);
            }
        }

        best_prob.unwrap()
    }
}

pub(crate) trait AltAlleleEmissionBuilder {
    type EmissionParams: stats::pairhmm::EmissionParameters + pairhmm::RefBaseEmission;

    fn build<'a>(
        &self,
        read_emission_params: &'a ReadEmission,
        ref_seq: &'a [u8],
    ) -> Self::EmissionParams;
}

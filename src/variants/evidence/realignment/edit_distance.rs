// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::Range;

use bio::alignment::AlignmentOperation;
use bio::pattern_matching::myers::{self, long};
use bio::stats::LogProb;

use crate::default_ref_base_emission;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::utils::homopolymers::HomopolymerIndelOperation;
use crate::variants::evidence::realignment::pairhmm::{RefBaseEmission, EDIT_BAND};

use super::pairhmm::{ReadVsAlleleEmission, VariantEmission};

#[derive(Debug, Clone, Copy, Derefable, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize)]
pub(crate) struct EditDistance(#[deref] pub(crate) u32);

impl EditDistance {
    pub(crate) fn update(&mut self, other_dist: &EditDistance) {
        self.0 += **other_dist;
    }
}

pub(crate) fn is_explainable_by_error_rates(
    n_subs: usize,
    n_ins: usize,
    n_del: usize,
    alignment_len: usize,
    alignment_properties: &AlignmentProperties,
    read_error_rate: LogProb,
) -> bool {
    let expected_count = |prob: LogProb| (alignment_len as f64) * prob.exp();
    let expected_n_subs = expected_count(read_error_rate);
    let expected_n_ins = expected_count(alignment_properties.gap_params.prob_insertion_artifact);
    let expected_n_del = expected_count(alignment_properties.gap_params.prob_deletion_artifact);

    n_subs as f64 <= expected_n_subs
        && n_ins as f64 <= expected_n_ins
        && n_del as f64 <= expected_n_del
}

#[derive(Debug, Clone, CopyGetters, Getters, PartialEq, Eq, Hash)]
pub(crate) struct EditOperationCounts {
    #[get_copy = "pub(crate)"]
    substitutions: usize,
    #[get_copy = "pub(crate)"]
    insertions: usize,
    #[get_copy = "pub(crate)"]
    deletions: usize,
    #[get_copy = "pub(crate)"]
    is_explainable_by_error_rates: bool,
    #[get_copy = "pub(crate)"]
    alignment_idx: usize,
    #[get = "pub(crate)"]
    ref_range: Range<usize>,
    #[get = "pub(crate)"]
    in_range_alignment: Vec<AlignmentOperation>,
}

impl EditOperationCounts {
    pub(crate) fn new(
        n_subs: usize,
        n_ins: usize,
        n_del: usize,
        start: usize,
        end: usize,
        alignment: Vec<AlignmentOperation>,
        alignment_properties: &AlignmentProperties,
        read_error_rate: LogProb,
        alignment_idx: usize,
    ) -> Self {
        assert!(start < end);
        let alignment_len = end - start;
        let explainable = is_explainable_by_error_rates(
            n_subs,
            n_ins,
            n_del,
            alignment_len,
            alignment_properties,
            read_error_rate,
        );

        EditOperationCounts {
            substitutions: n_subs,
            insertions: n_ins,
            deletions: n_del,
            ref_range: start..end,
            in_range_alignment: alignment,
            is_explainable_by_error_rates: explainable,
            alignment_idx,
        }
    }

    pub(crate) fn edit_distance(&self) -> EditDistance {
        EditDistance((self.substitutions + self.insertions + self.deletions) as u32)
    }
}

impl PartialOrd for EditOperationCounts {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.is_explainable_by_error_rates && !other.is_explainable_by_error_rates {
            return Some(Ordering::Less);
        } else if !self.is_explainable_by_error_rates && other.is_explainable_by_error_rates {
            return Some(Ordering::Greater);
        } else {
            Some(
                self.substitutions
                    .cmp(&other.substitutions)
                    .then(self.insertions.cmp(&other.insertions))
                    .then(self.deletions.cmp(&other.deletions)),
            )
        }
    }
}

impl Ord for EditOperationCounts {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

enum Myers {
    #[cfg(has_u128)]
    Short(myers::Myers<u128>),
    #[cfg(not(has_u128))]
    Short(myers::Myers<u64>),
    Long(long::Myers<u64>),
}

pub(crate) struct EditDistanceCalculation {
    myers: Myers,
    read_seq: Vec<u8>,
}

impl EditDistanceCalculation {
    pub(crate) fn max_pattern_len() -> usize {
        128
    }

    /// Create new instance.
    ///
    /// # Arguments
    /// * `read_seq` - read sequence in window (may not exceed 128 bases).
    pub(crate) fn new<P>(read_seq: P) -> Self
    where
        P: Iterator<Item = u8> + DoubleEndedIterator + ExactSizeIterator,
    {
        let l = read_seq.len();
        let read_seq = read_seq.collect();

        #[cfg(has_u128)]
        let num_bits = 128;
        #[cfg(not(has_u128))]
        let num_bits = 64;
        let myers = if l <= num_bits {
            Myers::Short(myers::Myers::new(&read_seq))
        } else {
            Myers::Long(long::Myers::new(&read_seq))
        };

        EditDistanceCalculation { myers, read_seq }
    }

    /// Returns a reasonable upper bound for the edit distance in order to band the pairHMM computation.
    /// We use the best edit distance and add 5.
    pub(crate) fn calc_best_hit(
        &mut self,
        emission_params: &ReadVsAlleleEmission,
        max_dist: Option<usize>,
        alignment_properties: &AlignmentProperties,
    ) -> Option<EditDistanceHit> {
        let ref_seq = || {
            (0..emission_params.len_x()).map(|i| emission_params.ref_base(i).to_ascii_uppercase())
        };
        let mut best_dist = usize::max_value();
        let mut positions = Vec::new();
        let max_dist = max_dist.unwrap_or(self.read_seq.len());

        let mut handle_match = |pos, dist: usize| match dist.cmp(&best_dist) {
            Ordering::Less => {
                positions.clear();
                positions.push(pos);
                best_dist = dist;
            }
            Ordering::Equal => {
                positions.push(pos);
            }
            Ordering::Greater => (),
        };

        let alignments: Vec<Alignment> = match &mut self.myers {
            Myers::Short(myers) => {
                let mut matches = myers.find_all_lazy(ref_seq(), max_dist as u8);
                for (pos, dist) in &mut matches {
                    handle_match(pos, dist as usize);
                }

                // collect alignments
                positions
                    .iter()
                    .cloned()
                    .map(|pos| {
                        let mut alignment = Alignment::default();
                        let (start, _) = matches.path_at(pos, &mut alignment.operations).unwrap();
                        alignment.start = start;
                        alignment
                    })
                    .collect()
            }
            Myers::Long(myers) => {
                let mut matches = myers.find_all_lazy(ref_seq(), max_dist);
                for (pos, dist) in &mut matches {
                    handle_match(pos, dist);
                }

                // collect alignments
                positions
                    .iter()
                    .cloned()
                    .map(|pos| {
                        let mut alignment = Alignment::default();
                        let (start, _) = matches.path_at(pos, &mut alignment.operations).unwrap();
                        // We find a pos relative to ref end, hence we have to project it to a position relative to
                        // the start.
                        alignment.start = start;
                        alignment
                    })
                    .collect()
            }
        };

        if positions.is_empty() {
            None
        } else {
            let start = alignments[0].start();
            // take the last (aka first because we are mapping backwards) position for an upper bound of the putative end
            let end = cmp::min(
                alignments.last().unwrap().start() + self.read_seq.len() + best_dist as usize,
                emission_params.len_x(),
            );

            // METHOD: obtain edit operations within the alt allele
            let edit_operation_counts = if emission_params.variant_ref_range().is_some() {
                Some(
                    alignments
                        .iter()
                        .enumerate()
                        .map(|(i, alignment)| {
                            let mut n_subst = 0;
                            let mut n_del = 0;
                            let mut n_ins = 0;
                            let mut pos = (emission_params.ref_offset() + alignment.start) as u64;
                            dbg!((alignment.start, emission_params.ref_offset()));
                            let pos_start = pos;
                            let mut in_range_alignment = Vec::new();
                            for op in alignment.operations() {
                                if emission_params.is_in_variant_ref_range(pos) {
                                    in_range_alignment.push(op.to_owned());
                                }
                                match op {
                                    AlignmentOperation::Subst => {
                                        if emission_params.is_in_variant_ref_range(pos) {
                                            n_subst += 1;
                                        }
                                        pos += 1;
                                    }
                                    AlignmentOperation::Del => {
                                        if emission_params.is_in_variant_ref_range(pos) {
                                            n_del += 1;
                                        }
                                        pos += 1;
                                    }
                                    AlignmentOperation::Ins => {
                                        if emission_params.is_in_variant_ref_range(pos) {
                                            n_ins += 1;
                                        }
                                    }
                                    AlignmentOperation::Match => {
                                        pos += 1;
                                    }
                                    _ => {
                                        unreachable!("bug: unexpected alignment operation");
                                    }
                                }
                            }

                            EditOperationCounts::new(
                                n_subst,
                                n_del,
                                n_ins,
                                pos_start as usize,
                                pos as usize,
                                in_range_alignment,
                                alignment_properties,
                                emission_params.read_emission().error_rate(),
                                i,
                            )
                        })
                        .min()
                        .unwrap(),
                )
            } else {
                None
            };

            // METHOD: obtain indel operations for homopolymer error model
            let homopolymer_indel_len = if emission_params.is_homopolymer_indel() {
                alignments
                    .iter()
                    .filter_map(|alignment| {
                        if let Some(operation) = HomopolymerIndelOperation::from_alignment(
                            &ref_seq().skip(alignment.start).collect::<Vec<_>>(),
                            &self.read_seq,
                            &alignment.operations,
                        ) {
                            if let Some(variant_ref_range) =
                                emission_params.variant_homopolymer_ref_range()
                            {
                                let ref_pos = (emission_params.ref_offset()
                                    + alignment.start
                                    + operation.text_pos())
                                    as u64;
                                // METHOD: check whether the operation is within the homopolymer variant range.
                                // In case of a deletion (operation.len() < 0) we also check whether the
                                // end of the deletion is within the variant ref range.
                                if variant_ref_range.contains(&(ref_pos))
                                    && (operation.len() > 0
                                        || variant_ref_range
                                            .contains(&(ref_pos + operation.len().abs() as u64)))
                                {
                                    Some(operation.len())
                                } else {
                                    None
                                }
                            } else {
                                None
                            }
                        } else {
                            None
                        }
                    })
                    .min_by_key(|indel_len| *indel_len)
            } else {
                None
            };

            Some(EditDistanceHit {
                start,
                end,
                dist: best_dist,
                alignments,
                homopolymer_indel_len,
                edit_operation_counts,
            })
        }
    }

    pub(crate) fn derive_allele_from_read<'a>(
        &self,
        emission_params: &'a ReadVsAlleleEmission,
        edit_distance_hit: &EditDistanceHit,
    ) -> Option<ReadVsAlleleEmission<'a>> {
        // Ensure that the given emission params are not shrunken, because the edit distance hit
        // is against the not shrunken parameters.
        assert_eq!(
            emission_params.ref_offset_orig(),
            emission_params.ref_offset()
        );
        assert_eq!(emission_params.ref_end_orig(), emission_params.ref_end());

        if edit_distance_hit
            .edit_operation_counts()
            .as_ref()
            .map_or(false, |opcounts| !opcounts.is_explainable_by_error_rates())
        {
            // METHOD: there are more errors than the error rates could explain
            // Therefore, we derive a novel allele from the read and consider this as an
            // alternative to the allele of interest. The read will map better to that,
            // such that it no longer counts as evidence for the allele of interest.
            // In case of ambiguity, this will also properly be captured by the comparison.

            // TODO in the future, this can be made more efficient by applying the patch transparently
            // with some helper data structure upon invocation of ref_base.

            let alignment = edit_distance_hit.some_alignment();

            // relative position in the emission snippet we have aligned against
            let mut pos_ref = alignment.start;
            let mut pos_read = 0; // semiglobal alignment
            let mut end_reduce = 0;

            let mut allele = Vec::new();
            let opcounts = edit_distance_hit.edit_operation_counts().as_ref().unwrap();
            // add part before the alignment
            allele.extend((0..alignment.start).map(|i| emission_params.ref_base(i)));
            dbg!((pos_ref, emission_params.ref_offset_orig()));

            for op in alignment.operations() {
                let is_in_range = emission_params
                    .is_in_variant_ref_range(pos_ref as u64 + emission_params.ref_offset() as u64);
                match op {
                    AlignmentOperation::Match => {
                        allele.push(emission_params.ref_base(pos_ref));
                        pos_ref += 1;
                        pos_read += 1;
                    }
                    AlignmentOperation::Subst => {
                        allele.push(if is_in_range {
                            self.read_seq[pos_read]
                        } else {
                            emission_params.ref_base(pos_ref)
                        });
                        pos_ref += 1;
                        pos_read += 1;
                    }
                    AlignmentOperation::Del => {
                        if is_in_range {
                            pos_ref += 1;
                            end_reduce += 1;
                        } else {
                            emission_params.ref_base(pos_ref);
                            pos_ref += 1;
                        }
                    }
                    AlignmentOperation::Ins => {
                        if is_in_range {
                            allele.push(self.read_seq[pos_read]);
                            pos_read += 1;
                        } else {
                            pos_read += 1;
                        }
                    }
                    _ => {
                        unreachable!("bug: unexpected alignment operation")
                    }
                }
            }
            // add the remaining sequence
            allele.extend(
                (pos_ref..emission_params.ref_end() - emission_params.ref_offset() + 1)
                    .map(|i| emission_params.ref_base(i)),
            );

            // Adjust ref_end by end_reduce, in order to ensure that deletions in var_range do not shorten the allele
            // without adjusting the ref_end.

            Some(ReadVsAlleleEmission::new(
                emission_params.read_emission(),
                Box::new(PatchedAlleleEmission {
                    ref_offset: emission_params.allele_emission().ref_offset(),
                    ref_end: emission_params.allele_emission().ref_end() - end_reduce,
                    patched_seq: allele,
                    ref_offset_override: None,
                    ref_end_override: None,
                }),
            ))
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Getters, CopyGetters, Default)]
pub(crate) struct Alignment {
    #[getset(get = "pub(crate)")]
    operations: Vec<AlignmentOperation>,
    #[getset(get_copy = "pub(crate)")]
    start: usize,
}

#[derive(Debug, Clone, Getters, CopyGetters)]
pub(crate) struct EditDistanceHit {
    #[getset(get_copy = "pub(crate)")]
    start: usize,
    #[getset(get_copy = "pub(crate)")]
    end: usize,
    #[getset(get_copy = "pub(crate)")]
    dist: usize,
    #[getset(get = "pub(crate)")]
    alignments: Vec<Alignment>,
    #[getset(get_copy = "pub(crate)")]
    homopolymer_indel_len: Option<i8>,
    #[getset(get = "pub(crate)")]
    edit_operation_counts: Option<EditOperationCounts>,
}

impl EditDistanceHit {
    pub(crate) fn dist_upper_bound(&self) -> usize {
        self.dist as usize + EDIT_BAND
    }

    /// Return the best alignment or any one if this cannot be decided.
    pub(crate) fn some_alignment(&self) -> &Alignment {
        if let Some(edit_operation_counts) = &self.edit_operation_counts {
            &self.alignments[edit_operation_counts.alignment_idx()]
        } else {
            &self.alignments[0]
        }
    }

    pub(crate) fn edit_distance(&self) -> Option<EditDistance> {
        self.edit_operation_counts()
            .as_ref()
            .map(|opcounts| opcounts.edit_distance())
    }
}

pub(crate) struct PatchedAlleleEmission {
    ref_offset: usize,
    ref_end: usize,
    ref_offset_override: Option<usize>,
    ref_end_override: Option<usize>,
    patched_seq: Vec<u8>,
}

impl RefBaseEmission for PatchedAlleleEmission {
    fn ref_base(&self, pos: usize) -> u8 {
        self.patched_seq[pos]
    }

    default_ref_base_emission!();

    fn variant_homopolymer_ref_range(&self) -> Option<Range<u64>> {
        unreachable!("bug: PatchedAlleleEmission should not be used as a normal alt allele");
    }

    fn variant_ref_range(&self) -> Option<Range<u64>> {
        unreachable!("bug: PatchedAlleleEmission should not be used as a normal alt allele");
    }

    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }
}

impl VariantEmission for PatchedAlleleEmission {
    fn is_homopolymer_indel(&self) -> bool {
        unreachable!("bug: PatchedAlleleEmission should not be used as a normal alt allele");
    }
}

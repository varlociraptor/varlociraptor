// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::cmp::Ordering;
use std::fmt::Debug;

use bio::alignment::AlignmentOperation;
use bio::pattern_matching::myers::{self, long};
use bio::stats::pairhmm;
use itertools::Itertools;

use crate::variants::evidence::realignment::pairhmm::{RefBaseEmission, EDIT_BAND};

enum Myers {
    Short(myers::Myers<u128>),
    Long(long::Myers<u64>),
}

pub(crate) struct EditDistanceCalculation {
    myers: Myers,
    read_seq_len: usize,
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

        let myers = if l <= 128 {
            Myers::Short(myers::Myers::new(read_seq))
        } else {
            Myers::Long(long::Myers::new(read_seq))
        };

        EditDistanceCalculation {
            myers,
            read_seq_len: l,
        }
    }

    /// Returns a reasonable upper bound for the edit distance in order to band the pairHMM computation.
    /// We use the best edit distance and add 5.
    pub(crate) fn calc_best_hit<E: pairhmm::EmissionParameters + RefBaseEmission>(
        &mut self,
        emission_params: &E,
        max_dist: Option<usize>,
    ) -> Option<EditDistanceHit> {
        let ref_seq =
            (0..emission_params.len_x()).map(|i| emission_params.ref_base(i).to_ascii_uppercase());
        let mut best_dist = usize::max_value();
        let mut positions = Vec::new();
        let alignments: Vec<Alignment>;
        let max_dist = max_dist.unwrap_or(self.read_seq_len);

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

        match &mut self.myers {
            Myers::Short(myers) => {
                let mut matches = myers.find_all_lazy(ref_seq, max_dist as u8);
                for (pos, dist) in &mut matches {
                    handle_match(pos, dist as usize);
                }

                // collect alignments
                alignments = positions
                    .iter()
                    .cloned()
                    .map(|pos| {
                        let mut alignment = Alignment::default();
                        let (start, _) = matches.path_at(pos, &mut alignment.operations).unwrap();
                        alignment.start = start;
                        alignment
                    })
                    .collect();
            }
            Myers::Long(myers) => {
                let mut matches = myers.find_all_lazy(ref_seq, max_dist);
                for (pos, dist) in &mut matches {
                    handle_match(pos, dist);
                }

                // collect alignments
                alignments = positions
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
                    .collect();
            }
        }

        if positions.is_empty() {
            None
        } else {
            let start = alignments[0].start();
            // take the last (aka first because we are mapping backwards) position for an upper bound of the putative end
            let end = cmp::min(
                alignments.last().unwrap().start() + self.read_seq_len + best_dist as usize,
                emission_params.len_x(),
            );

            // METHOD: obtain indel operations for divindel bias
            let best_indel_operations = alignments
                .iter()
                .map(|alignment| {
                    let mut ref_pos = emission_params.ref_offset() + alignment.start;

                    // METHOD: group operations by indel or not, check whether they overlap the variant,
                    // and report additional indels that remain in the variant range.
                    // These are indicative of divindel bias (i.e. some wild disagreeing indel operations cause a caller
                    // to interpet them as an indel variant, but in reality it is e.g. a PCR homopolymer error.)
                    alignment
                        .operations()
                        .iter()
                        .group_by(|op| {
                            matches!(op, AlignmentOperation::Del | AlignmentOperation::Ins)
                        })
                        .into_iter()
                        .filter_map(|(is_indel, group)| {
                            let group: Vec<_> = group.collect();

                            let mut group_ref_len = 0;
                            for op in &group {
                                // update ref_pos
                                match op {
                                    AlignmentOperation::Match
                                    | AlignmentOperation::Subst
                                    | AlignmentOperation::Del => group_ref_len += 1,
                                    AlignmentOperation::Ins => (),
                                    _ => unreachable!(),
                                }
                            }

                            let ret = if let Some(variant_ref_range) =
                                emission_params.variant_ref_range()
                            {
                                if is_indel
                                    && (variant_ref_range.contains(&ref_pos)
                                        || variant_ref_range.contains(&(ref_pos + group_ref_len))
                                        || (variant_ref_range.start >= ref_pos
                                            && variant_ref_range.end <= ref_pos + group_ref_len))
                                {
                                    // METHOD: indel ops that overlap with the variant interval are recorded here.
                                    Some(group)
                                } else {
                                    None
                                }
                            } else {
                                None
                            };

                            ref_pos += group_ref_len;

                            ret
                        })
                        .flatten()
                        .cloned()
                        .collect_vec()
                })
                .min_by_key(|indels| indels.len())
                .unwrap();

            Some(EditDistanceHit {
                start,
                end,
                dist: best_dist,
                alignments,
                best_indel_operations,
            })
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
    #[getset(get = "pub(crate)")]
    best_indel_operations: Vec<AlignmentOperation>,
}

impl EditDistanceHit {
    pub(crate) fn dist_upper_bound(&self) -> usize {
        self.dist as usize + EDIT_BAND
    }
}

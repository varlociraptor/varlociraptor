// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::cmp::Ordering;
use std::fmt::Debug;

use bio::alignment::AlignmentOperation;
use bio::pattern_matching::myers::{self, long};

use crate::utils::homopolymers::HomopolymerIndelOperation;
use crate::variants::evidence::realignment::pairhmm::{RefBaseEmission, EDIT_BAND};

use super::pairhmm::{ReadVsAlleleEmission, VariantEmission};

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
    ) -> Option<EditDistanceHit> {
        let ref_seq = || {
            (0..emission_params.len_x()).map(|i| emission_params.ref_base(i).to_ascii_uppercase())
        };
        let mut best_dist = usize::max_value();
        let mut positions = Vec::new();
        let alignments: Vec<Alignment>;
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

        match &mut self.myers {
            Myers::Short(myers) => {
                let mut matches = myers.find_all_lazy(ref_seq(), max_dist as u8);
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
                let mut matches = myers.find_all_lazy(ref_seq(), max_dist);
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
                alignments.last().unwrap().start() + self.read_seq.len() + best_dist as usize,
                emission_params.len_x(),
            );

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
    #[getset(get_copy = "pub(crate)")]
    homopolymer_indel_len: Option<i8>,
}

impl EditDistanceHit {
    pub(crate) fn dist_upper_bound(&self) -> usize {
        self.dist as usize + EDIT_BAND
    }
}

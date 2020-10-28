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
            Some(EditDistanceHit {
                start,
                end,
                dist: best_dist,
                alignments,
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
}

impl EditDistanceHit {
    pub(crate) fn dist_upper_bound(&self) -> usize {
        self.dist as usize + EDIT_BAND
    }
}

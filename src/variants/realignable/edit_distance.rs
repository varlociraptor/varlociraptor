use std::cmp;
use std::fmt::Debug;

use anyhow::Result;
use bio::pattern_matching::myers::Myers;
use bio::stats::pairhmm;
use bio::stats::{LogProb, PHREDProb, Prob};
use bio_types::strand::Strand;
use rust_htslib::bam;
use rust_htslib::bam::record::CigarStringView;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::model::Variant;
use crate::variants::realignable::pairhmm::{RefBaseEmission, EDIT_BAND};

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

#[derive(Debug, Clone, CopyGetters)]
#[getset(get_copy = "pub")]
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

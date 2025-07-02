use bio_types::sequence::SequenceReadPairOrientation;

use crate::variants::evidence::observations::depth_observation::ProcessedDepthObservation;

use super::{depth_observation::DepthObservation, read_observation::ProcessedReadObservation};

#[derive(Debug, Getters, MutGetters, Default)]
#[getset(get = "pub", get_mut = "pub(crate)")]
pub struct Pileup {
    read_observations: Vec<ProcessedReadObservation>,
    depth_observations: Vec<ProcessedDepthObservation>,
    n_filtered_out_observations: usize,
}

impl Pileup {
    pub(crate) fn new(
        read_observations: Vec<ProcessedReadObservation>,
        depth_observations: Vec<ProcessedDepthObservation>,
    ) -> Self {
        Self {
            read_observations,
            depth_observations,
            n_filtered_out_observations: 0,
        }
    }

    /// Remove all non-standard alignments from pileup (softclipped observations, non-standard read orientations).
    pub(crate) fn remove_nonstandard_alignments(
        &mut self,
        omit_read_orientation_bias: bool,
    ) -> bool {
        // METHOD: this can be helpful to get cleaner SNV and MNV calls. Support for those should be
        // solely driven by standard alignments, that are in expected orientation.
        // Otherwise called SNVs can be artifacts of near SVs.
        let n_orig = self.read_observations.len();
        self.read_observations.retain(|obs| {
            omit_read_orientation_bias
                || (obs.read_orientation == SequenceReadPairOrientation::F1R2
                    || obs.read_orientation == SequenceReadPairOrientation::F2R1
                    || obs.read_orientation == SequenceReadPairOrientation::None)
        });
        self.n_filtered_out_observations += n_orig - self.read_observations.len();

        self.n_filtered_out_observations > 0
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.read_observations.is_empty() && self.depth_observations.is_empty()
    }
}

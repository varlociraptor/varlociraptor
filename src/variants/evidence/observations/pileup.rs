use bio_types::sequence::SequenceReadPairOrientation;

use super::{depth_observation::DepthObservation, read_observation::ProcessedReadObservation};

#[derive(new, Debug, Getters, MutGetters, Default)]
#[getset(get = "pub(crate)", get_mut = "pub(crate)")]
pub(crate) struct Pileup {
    read_observations: Vec<ProcessedReadObservation>,
    depth_observations: Vec<DepthObservation>,
}

impl Pileup {
    /// Remove all non-standard alignments from pileup (softclipped observations, non-standard read orientations).
    pub(crate) fn remove_nonstandard_alignments(&mut self, omit_read_orientation_bias: bool) {
        // METHOD: this can be helpful to get cleaner SNV and MNV calls. Support for those should be
        // solely driven by standard alignments, that are in expected orientation.
        // Otherwise called SNVs can be artifacts of near SVs.
        self.read_observations.retain(|obs| {
            omit_read_orientation_bias
                || (obs.read_orientation == SequenceReadPairOrientation::F1R2
                    || obs.read_orientation == SequenceReadPairOrientation::F2R1
                    || obs.read_orientation == SequenceReadPairOrientation::None)
        });
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.read_observations.is_empty() && self.depth_observations.is_empty()
    }
}

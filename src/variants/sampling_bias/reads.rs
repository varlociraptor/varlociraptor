// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;

use bio::stats::LogProb;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::sampling_bias::SamplingBias;
use crate::variants::types::Variant;

pub(crate) trait ReadSamplingBias: Variant + SamplingBias {
    /// Probability to sample read from alt allele given the average feasible positions observed
    /// from a subsample of the mapped reads.
    ///
    /// The key idea is calculate the probability as number of valid placements (considering the
    /// max softclip allowed by the mapper) over all possible placements.
    fn prob_sample_alt_read(
        &self,
        read_len: u64,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        if let Some(feasible) = self.feasible_bases(read_len, alignment_properties) {
            let prob = {
                let n_alt = self
                    .enclosable_len()
                    .map_or(read_len, |len| cmp::min(len, read_len));
                let n_alt_valid = cmp::min(n_alt, feasible);

                LogProb((n_alt_valid as f64).ln() - (n_alt as f64).ln())
            };
            assert!(prob.is_valid());

            dbg!(prob);

            prob
        } else {
            LogProb::ln_one()
        }
    }
}

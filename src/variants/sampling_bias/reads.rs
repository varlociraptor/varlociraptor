use std::cmp;

use bio::stats::LogProb;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::sampling_bias::SamplingBias;
use crate::variants::Variant;

pub trait ReadSamplingBias<'a>: Variant<'a> + SamplingBias<'a> {
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
        let feasible = self.feasible_bases(read_len, alignment_properties);

        let prob = {
            let n_alt = cmp::min(self.len(), read_len);
            let n_alt_valid = cmp::min(n_alt, feasible);

            LogProb((n_alt_valid as f64).ln() - (n_alt as f64).ln())
        };
        assert!(prob.is_valid());

        prob
    }
}

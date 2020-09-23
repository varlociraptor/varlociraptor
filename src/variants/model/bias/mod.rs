use std::cmp::Ordering;

use bio::stats::probs::LogProb;
use itertools::Itertools;
use strum::IntoEnumIterator;

use crate::variants::evidence::observation::Observation;

pub(crate) mod read_orientation_bias;
pub(crate) mod strand_bias;

pub(crate) use read_orientation_bias::ReadOrientationBias;
pub(crate) use strand_bias::StrandBias;

pub(crate) trait Bias {
    fn prob(&self, observation: &Observation) -> LogProb;

    fn prob_any(&self) -> LogProb;

    fn is_artifact(&self) -> bool;
}

#[derive(Builder, CopyGetters, Getters, Debug, Clone)]
#[builder(build_fn(name = "build_inner"))]
pub(crate) struct Biases {
    #[getset(get = "pub(crate)")]
    strand_bias: StrandBias,
    #[getset(get = "pub(crate)")]
    read_orientation_bias: ReadOrientationBias,
    #[getset(get_copy = "pub(crate)")]
    #[builder(private)]
    pub(crate) prob_any: LogProb,
}

impl PartialEq for Biases {
    fn eq(&self, other: &Self) -> bool {
        self.strand_bias == other.strand_bias
            && self.read_orientation_bias == other.read_orientation_bias
    }
}

impl Eq for Biases {}

impl Ord for Biases {
    fn cmp(&self, other: &Self) -> Ordering {
        self.strand_bias
            .cmp(&other.strand_bias)
            .then(self.read_orientation_bias.cmp(&other.read_orientation_bias))
    }
}

impl PartialOrd for Biases {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl BiasesBuilder {
    pub(crate) fn build(&mut self) -> Result<Biases, String> {
        let prob_any = self
            .strand_bias
            .expect("bug: strand_bias must be set before building")
            .prob_any()
            + self
                .read_orientation_bias
                .expect("bug: read_orientation_bias must be set before building")
                .prob_any();

        self.prob_any(prob_any);
        self.build_inner()
    }
}

impl Biases {
    pub(crate) fn all_artifact_combinations(
        consider_read_orientation_bias: bool,
    ) -> Box<dyn Iterator<Item = Self>> {
        if consider_read_orientation_bias {
            Box::new(
                StrandBias::iter()
                    .cartesian_product(ReadOrientationBias::iter())
                    .map(|(strand_bias, read_orientation_bias)| {
                        BiasesBuilder::default()
                            .strand_bias(strand_bias)
                            .read_orientation_bias(read_orientation_bias)
                            .build()
                            .unwrap()
                    }),
            )
        } else {
            Box::new(StrandBias::iter().map(|strand_bias| {
                BiasesBuilder::default()
                    .strand_bias(strand_bias)
                    .read_orientation_bias(ReadOrientationBias::None)
                    .build()
                    .unwrap()
            }))
        }
    }

    pub(crate) fn none() -> Self {
        BiasesBuilder::default()
            .strand_bias(StrandBias::None)
            .read_orientation_bias(ReadOrientationBias::None)
            .build()
            .unwrap()
    }

    pub(crate) fn prob(&self, observation: &Observation) -> LogProb {
        self.strand_bias.prob(observation) + self.read_orientation_bias.prob(observation)
    }

    pub(crate) fn is_artifact(&self) -> bool {
        self.strand_bias.is_artifact() || self.read_orientation_bias.is_artifact()
    }
}

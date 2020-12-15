use std::cmp::Ordering;

use bio::stats::probs::LogProb;
use itertools::Itertools;
use strum::IntoEnumIterator;

use crate::variants::evidence::observation::{Observation, ReadPosition};

pub(crate) mod read_orientation_bias;
pub(crate) mod read_position_bias;
pub(crate) mod strand_bias;

pub(crate) use read_orientation_bias::ReadOrientationBias;
pub(crate) use read_position_bias::ReadPositionBias;
pub(crate) use strand_bias::StrandBias;

pub(crate) trait Bias {
    fn prob(&self, observation: &Observation<ReadPosition>) -> LogProb;

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
    #[getset(get = "pub(crate)")]
    read_position_bias: ReadPositionBias,
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
        consider_strand_bias: bool,
        consider_read_position_bias: bool,
    ) -> Box<dyn Iterator<Item = Self>> {
        if !consider_strand_bias && !consider_read_orientation_bias && !consider_read_position_bias
        {
            return Box::new(std::iter::empty());
        }

        let strand_biases = if consider_strand_bias {
            StrandBias::iter().collect_vec()
        } else {
            vec![StrandBias::None]
        };
        let read_position_biases = if consider_read_position_bias {
            ReadPositionBias::iter().collect_vec()
        } else {
            vec![ReadPositionBias::None]
        };
        let read_orientation_biases = if consider_read_orientation_bias {
            ReadOrientationBias::iter().collect_vec()
        } else {
            vec![ReadOrientationBias::None]
        };

        Box::new(
            strand_biases
                .into_iter()
                .cartesian_product(read_orientation_biases.into_iter())
                .cartesian_product(read_position_biases.into_iter())
                .filter_map(|((sb, rob), rpb)| {
                    if sb.is_artifact() || rob.is_artifact() || rpb.is_artifact() {
                        Some(
                            BiasesBuilder::default()
                                .strand_bias(sb)
                                .read_orientation_bias(rob)
                                .read_position_bias(rpb)
                                .build()
                                .unwrap(),
                        )
                    } else {
                        None
                    }
                }),
        )
    }

    pub(crate) fn none() -> Self {
        BiasesBuilder::default()
            .strand_bias(StrandBias::None)
            .read_orientation_bias(ReadOrientationBias::None)
            .read_position_bias(ReadPositionBias::None)
            .build()
            .unwrap()
    }

    pub(crate) fn prob(&self, observation: &Observation<ReadPosition>) -> LogProb {
        self.strand_bias.prob(observation) + self.read_orientation_bias.prob(observation)
    }

    pub(crate) fn is_artifact(&self) -> bool {
        self.strand_bias.is_artifact() || self.read_orientation_bias.is_artifact()
    }
}

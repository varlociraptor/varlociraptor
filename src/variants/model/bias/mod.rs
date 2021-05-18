

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

    fn prob_any(&self, observation: &Observation<ReadPosition>) -> LogProb;

    fn is_artifact(&self) -> bool;

    fn is_possible(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        pileups.iter().any(|pileup| {
            pileup
                .iter()
                .any(|observation| self.prob(observation) != LogProb::ln_zero())
        })
    }

    #[allow(unused_variables)]
    fn is_informative(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        true
    }
}

#[derive(Builder, CopyGetters, Getters, Debug, Clone, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub(crate) struct Biases {
    #[getset(get = "pub(crate)")]
    strand_bias: StrandBias,
    #[getset(get = "pub(crate)")]
    read_orientation_bias: ReadOrientationBias,
    #[getset(get = "pub(crate)")]
    read_position_bias: ReadPositionBias,
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
                    match (sb.is_artifact(), rob.is_artifact(), rpb.is_artifact()) {
                        (true, false, false) | (false, true, false) | (false, false, true) => Some(
                            BiasesBuilder::default()
                                .strand_bias(sb)
                                .read_orientation_bias(rob)
                                .read_position_bias(rpb)
                                .build()
                                .unwrap(),
                        ),
                        _ => None,
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

    pub(crate) fn is_possible(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        self.strand_bias.is_possible(pileups)
            && self.read_orientation_bias.is_possible(pileups)
            && self.read_position_bias.is_possible(pileups)
    }

    pub(crate) fn is_informative(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        self.strand_bias.is_informative(pileups)
            && self.read_orientation_bias.is_informative(pileups)
            && self.read_position_bias.is_informative(pileups)
    }

    pub(crate) fn prob(&self, observation: &Observation<ReadPosition>) -> LogProb {
        //dbg!(self.strand_bias.prob(observation), self.read_orientation_bias.prob(observation), self.read_position_bias.prob(observation));
        self.strand_bias.prob(observation)
            + self.read_orientation_bias.prob(observation)
            + self.read_position_bias.prob(observation)
    }

    pub(crate) fn prob_any(&self, observation: &Observation<ReadPosition>) -> LogProb {
        self.strand_bias.prob_any(observation)
            + self.read_orientation_bias.prob_any(observation)
            + self.read_position_bias.prob_any(observation)
    }

    pub(crate) fn is_artifact(&self) -> bool {
        self.strand_bias.is_artifact()
            || self.read_orientation_bias.is_artifact()
            || self.read_position_bias.is_artifact()
    }
}

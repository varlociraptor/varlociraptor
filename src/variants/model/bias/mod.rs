use std::cmp;

use bio::stats::bayesian::bayes_factors::{evidence::KassRaftery, BayesFactor};
use bio::stats::probs::LogProb;
use itertools::Itertools;
use strum::IntoEnumIterator;

use crate::utils::PROB_095;
use crate::variants::evidence::observation::{IndelOperations, Observation, ReadPosition};

pub(crate) mod divindel_bias;
pub(crate) mod read_orientation_bias;
pub(crate) mod read_position_bias;
pub(crate) mod softclip_bias;
pub(crate) mod strand_bias;

pub(crate) use divindel_bias::DivIndelBias;
pub(crate) use read_orientation_bias::ReadOrientationBias;
pub(crate) use read_position_bias::ReadPositionBias;
pub(crate) use softclip_bias::SoftclipBias;
pub(crate) use strand_bias::StrandBias;

pub(crate) trait Bias: Default + cmp::PartialEq {
    fn prob(&self, observation: &Observation<ReadPosition, IndelOperations>) -> LogProb;

    fn prob_any(&self, observation: &Observation<ReadPosition, IndelOperations>) -> LogProb;

    fn is_artifact(&self) -> bool;

    fn is_possible(&self, pileups: &[Vec<Observation<ReadPosition, IndelOperations>>]) -> bool {
        pileups.iter().any(|pileup| {
            pileup
                .iter()
                .any(|observation| self.prob(observation) != LogProb::ln_zero())
        })
    }

    fn is_informative(&self, _pileups: &[Vec<Observation<ReadPosition, IndelOperations>>]) -> bool {
        true
    }

    fn is_likely(&self, pileups: &[Vec<Observation<ReadPosition, IndelOperations>>]) -> bool {
        if *self == Self::default() {
            true
        } else {
            pileups.iter().any(|pileup| {
                let strong_all = pileup.iter().filter(&Self::is_strong_obs).count();
                if strong_all >= 10 {
                    let strong_bias_evidence = pileup
                        .iter()
                        .filter(|obs| {
                            Self::is_strong_obs(obs) && self.prob(obs) == LogProb::ln_one()
                        })
                        .count();
                    // METHOD: there is bias evidence if we have at least two third of the strong observations supporting the bias
                    let ratio = strong_bias_evidence as f64 / strong_all as f64;
                    ratio >= 0.66666
                } else {
                    // METHOD: not enough reads, rather consider all biases to be sure
                    true
                }
            })
        }
    }

    /// Learn parameters needed for estimation on current pileup.
    fn learn_parameters(&mut self, _pileups: &[Vec<Observation<ReadPosition, IndelOperations>>]) {
        // METHOD: by default, there is nothing to learn, however, a bias can use this to
        // infer some parameters over which we would otherwise need to integrate (which would hamper
        // performance too much).
    }

    fn is_strong_obs(obs: &&Observation<ReadPosition, IndelOperations>) -> bool {
        obs.prob_mapping() >= *PROB_095
            && BayesFactor::new(obs.prob_alt, obs.prob_ref).evidence_kass_raftery()
                >= KassRaftery::Strong
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
    #[getset(get = "pub(crate)")]
    softclip_bias: SoftclipBias,
    #[getset(get = "pub(crate)")]
    divindel_bias: DivIndelBias,
}

impl Biases {
    pub(crate) fn all_artifact_combinations(
        consider_read_orientation_bias: bool,
        consider_strand_bias: bool,
        consider_read_position_bias: bool,
        consider_softclip_bias: bool,
        consider_divindel_bias: bool,
    ) -> Box<dyn Iterator<Item = Self>> {
        if !consider_strand_bias
            && !consider_read_orientation_bias
            && !consider_read_position_bias
            && !consider_softclip_bias
            && !consider_divindel_bias
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
        let softclip_biases = if consider_softclip_bias {
            SoftclipBias::iter().collect_vec()
        } else {
            vec![SoftclipBias::None]
        };
        let divindel_biases = if consider_divindel_bias {
            DivIndelBias::values()
        } else {
            vec![DivIndelBias::None]
        };

        Box::new(
            strand_biases
                .into_iter()
                .cartesian_product(read_orientation_biases.into_iter())
                .cartesian_product(read_position_biases.into_iter())
                .cartesian_product(softclip_biases.into_iter())
                .cartesian_product(divindel_biases.into_iter())
                .filter_map(|((((sb, rob), rpb), scb), dib)| {
                    if [
                        sb.is_artifact(),
                        rob.is_artifact(),
                        rpb.is_artifact(),
                        scb.is_artifact(),
                        dib.is_artifact(),
                    ]
                    .iter()
                    .map(|artifact| if *artifact { 1 } else { 0 })
                    .sum::<usize>()
                        == 1
                    {
                        Some(
                            BiasesBuilder::default()
                                .strand_bias(sb)
                                .read_orientation_bias(rob)
                                .read_position_bias(rpb)
                                .softclip_bias(scb)
                                .divindel_bias(dib)
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
            .softclip_bias(SoftclipBias::None)
            .divindel_bias(DivIndelBias::None)
            .build()
            .unwrap()
    }

    pub(crate) fn is_possible(
        &self,
        pileups: &[Vec<Observation<ReadPosition, IndelOperations>>],
    ) -> bool {
        self.strand_bias.is_possible(pileups)
            && self.read_orientation_bias.is_possible(pileups)
            && self.read_position_bias.is_possible(pileups)
            && self.softclip_bias.is_possible(pileups)
            && self.divindel_bias.is_possible(pileups)
    }

    pub(crate) fn is_informative(
        &self,
        pileups: &[Vec<Observation<ReadPosition, IndelOperations>>],
    ) -> bool {
        self.strand_bias.is_informative(pileups)
            && self.read_orientation_bias.is_informative(pileups)
            && self.read_position_bias.is_informative(pileups)
            && self.softclip_bias.is_informative(pileups)
            && self.divindel_bias.is_informative(pileups)
    }

    pub(crate) fn is_likely(
        &self,
        pileups: &[Vec<Observation<ReadPosition, IndelOperations>>],
    ) -> bool {
        self.strand_bias.is_likely(pileups)
            && self.read_orientation_bias.is_likely(pileups)
            && self.read_position_bias.is_likely(pileups)
            && self.softclip_bias.is_likely(pileups)
            && self.divindel_bias.is_likely(pileups)
    }

    pub(crate) fn prob(&self, observation: &Observation<ReadPosition, IndelOperations>) -> LogProb {
        //dbg!(self.strand_bias.prob(observation), self.read_orientation_bias.prob(observation), self.read_position_bias.prob(observation));
        self.strand_bias.prob(observation)
            + self.read_orientation_bias.prob(observation)
            + self.read_position_bias.prob(observation)
            + self.softclip_bias.prob(observation)
            + self.divindel_bias.prob(observation)
    }

    pub(crate) fn prob_any(
        &self,
        observation: &Observation<ReadPosition, IndelOperations>,
    ) -> LogProb {
        self.strand_bias.prob_any(observation)
            + self.read_orientation_bias.prob_any(observation)
            + self.read_position_bias.prob_any(observation)
            + self.softclip_bias.prob_any(observation)
            + self.divindel_bias.prob_any(observation)
    }

    pub(crate) fn is_artifact(&self) -> bool {
        self.strand_bias.is_artifact()
            || self.read_orientation_bias.is_artifact()
            || self.read_position_bias.is_artifact()
            || self.softclip_bias.is_artifact()
            || self.divindel_bias.is_artifact()
    }

    pub(crate) fn learn_parameters(
        &mut self,
        pileups: &[Vec<Observation<ReadPosition, IndelOperations>>],
    ) {
        self.divindel_bias.learn_parameters(pileups)
    }
}

use bio::stats::probs::LogProb;

use crate::variants::evidence::observations::read_observation::{
    ProcessedReadObservation, ReadPosition,
};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash)]
pub(crate) enum ReadPositionBias {
    None,
    Some,
}

impl Default for ReadPositionBias {
    fn default() -> Self {
        ReadPositionBias::None
    }
}

impl Bias for ReadPositionBias {
    fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
        match (self, observation.read_position) {
            (ReadPositionBias::None, _) => observation.prob_hit_base, // normal
            (ReadPositionBias::Some, ReadPosition::Major) => LogProb::ln_one(), // bias
            (ReadPositionBias::Some, ReadPosition::Some) => LogProb::ln_zero(), // no bias
        }
    }

    fn prob_any(&self, observation: &ProcessedReadObservation) -> LogProb {
        observation.prob_hit_base
    }

    fn is_artifact(&self) -> bool {
        *self != ReadPositionBias::None
    }
}

// use bio::stats::probs::LogProb;
// use bio::stats::Prob;
// use itertools::Itertools;
// use ordered_float::NotNan;

// use crate::variants::evidence::observations::pileup::Pileup;
// use crate::variants::evidence::observations::read_observation::{
//     ProcessedReadObservation, ReadPosition,
// };
// use crate::variants::model::bias::Bias;

// #[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash)]
// pub(crate) enum ReadPositionBias {
//     None { major_rate: Option<NotNan<f64>> },
//     Some,
// }

// impl Default for ReadPositionBias {
//     fn default() -> Self {
//         ReadPositionBias::None { major_rate: None }
//     }
// }

// impl Bias for ReadPositionBias {
//     fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
//         match (self, observation.read_position) {
//             (
//                 ReadPositionBias::None {
//                     major_rate: Some(rate),
//                 },
//                 ReadPosition::Major,
//             ) => LogProb::from(Prob(**rate)),
//             (
//                 ReadPositionBias::None {
//                     major_rate: Some(rate),
//                 },
//                 ReadPosition::Some,
//             ) => LogProb::from(Prob(1.0 - **rate)),
//             (ReadPositionBias::None { major_rate: None }, ReadPosition::Major) => {
//                 observation.prob_hit_base
//             }
//             (ReadPositionBias::None { major_rate: None }, ReadPosition::Some) => {
//                 observation.prob_hit_base.ln_one_minus_exp()
//             }
//             (ReadPositionBias::Some, ReadPosition::Major) => LogProb::ln_one(), // bias
//             (ReadPositionBias::Some, ReadPosition::Some) => LogProb::ln_zero(), // no bias
//         }
//     }

//     fn prob_any(&self, observation: &ProcessedReadObservation) -> LogProb {
//         LogProb::ln_one()
//     }

//     fn is_artifact(&self) -> bool {
//         !matches!(self, ReadPositionBias::None { .. })
//     }

//     fn is_informative(&self, pileups: &[Pileup]) -> bool {
//         // METHOD: if all reads overlap the variant at the major pos,
//         // we cannot estimate a read position bias and None is the only informative one.
//         !self.is_artifact() || Self::estimate_major_rate(pileups).is_some()
//     }

//     fn learn_parameters(&mut self, pileups: &[Pileup]) {
//         if let ReadPositionBias::None { ref mut major_rate } = self {
//             *major_rate = Self::estimate_major_rate(pileups);
//         }
//     }
// }

// impl ReadPositionBias {
//     fn estimate_major_rate(pileups: &[Pileup]) -> Option<NotNan<f64>> {
//         let strong_all = LogProb::ln_sum_exp(
//             &pileups
//                 .iter()
//                 .flat_map(|pileup| {
//                     pileup.read_observations().iter().filter_map(|obs| {
//                         if obs.is_strong_ref_support() {
//                             Some(obs.prob_mapping())
//                         } else {
//                             None
//                         }
//                     })
//                 })
//                 .collect_vec(),
//         )
//         .exp();
//         let strong_major = LogProb::ln_sum_exp(
//             &pileups
//                 .iter()
//                 .flat_map(|pileup| {
//                     pileup.read_observations().iter().filter_map(|obs| {
//                         if obs.is_strong_ref_support() && obs.read_position == ReadPosition::Major {
//                             Some(obs.prob_mapping())
//                         } else {
//                             None
//                         }
//                     })
//                 })
//                 .collect_vec(),
//         )
//         .exp();
//         let any_major = pileups.iter().any(|pileup| {
//             pileup
//                 .read_observations()
//                 .iter()
//                 .any(|obs| obs.read_position == ReadPosition::Major)
//         });

//         if strong_all > 2.0 {
//             let major_fraction = strong_major / strong_all;
//             if any_major && major_fraction < 1.0 {
//                 // METHOD: if there is any read with the major read position in either
//                 // the ref or the alt supporting strong evidences and not all the reads
//                 // supporting ref do that at the major position, we report a fraction
//                 // and consider the read position bias.
//                 return Some(NotNan::new(major_fraction).unwrap());
//             }
//         }
//         None
//     }
// }

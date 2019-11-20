// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Deref, Range};

use bio::stats::LogProb;
use itertools::Itertools;
use ordered_float::NotNan;
use strum_macros::{EnumIter, EnumString, IntoStaticStr};

use crate::grammar;

pub mod evidence;
pub mod likelihood;
pub mod modes;
pub mod sample;

#[derive(Debug, Clone)]
pub struct Contamination {
    pub by: usize,
    pub fraction: f64,
}

#[derive(Ord, Eq, PartialOrd, PartialEq, Clone, Debug)]
pub struct Event {
    pub name: String,
    pub vafs: grammar::VAFTree,
    pub strand_bias: StrandBias,
}

impl Event {
    pub fn is_artifact(&self) -> bool {
        self.strand_bias != StrandBias::None
    }
}

pub type AlleleFreq = NotNan<f64>;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord)]
pub enum StrandBias {
    None,
    Forward,
    Reverse,
}

impl Default for StrandBias {
    fn default() -> Self {
        StrandBias::None
    }
}

impl StrandBias {
    pub fn is_some(&self) -> bool {
        if let StrandBias::None = self {
            false
        } else {
            true
        }
    }

    pub fn forward_rate(&self) -> LogProb {
        match self {
            StrandBias::None => LogProb(0.5_f64.ln()),
            StrandBias::Forward => LogProb::ln_one(),
            StrandBias::Reverse => LogProb::ln_zero(),
        }
    }

    pub fn reverse_rate(&self) -> LogProb {
        match self {
            StrandBias::None => LogProb(0.5_f64.ln()),
            StrandBias::Forward => LogProb::ln_zero(),
            StrandBias::Reverse => LogProb::ln_one(),
        }
    }
}

#[allow(non_snake_case)]
pub fn AlleleFreq(af: f64) -> AlleleFreq {
    NotNan::new(af).unwrap()
}

pub trait AlleleFreqs: Debug {
    fn is_absent(&self) -> bool;
}
impl AlleleFreqs for DiscreteAlleleFreqs {
    fn is_absent(&self) -> bool {
        self.inner.len() == 1 && self.inner[0] == AlleleFreq(0.0)
    }
}
impl AlleleFreqs for ContinuousAlleleFreqs {
    fn is_absent(&self) -> bool {
        self.is_singleton() && self.start == AlleleFreq(0.0)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct DiscreteAlleleFreqs {
    inner: Vec<AlleleFreq>,
}

impl DiscreteAlleleFreqs {
    /// Create spectrum of discrete allele frequencies with given values.
    pub fn new(spectrum: Vec<AlleleFreq>) -> Self {
        DiscreteAlleleFreqs { inner: spectrum }
    }

    /// Return spectrum of all feasible allele frequencies for a given ploidy and maximum
    /// number of amplification.
    ///
    /// In principle, a ploidy of, e.g., 2 allows allele frequencies 0, 0.5, 1.0. However,
    /// if a locus is amplified e.g. 1 time, allele frequencies can be effectively
    /// 0, 0.25, 0.5, 0.75, 1.0, because all reads are projected to the original location in the
    /// reference genome, and it is unclear whether the first, the second, or both loci contain
    /// the variant.
    ///
    /// # Arguments
    /// * ploidy - the assumed overall ploidy
    /// * max_amplification - the maximum amplification factor (1 means no amplification, 2 means
    ///   at most one duplicate, ...).
    fn _feasible(ploidy: u32, max_amplification: u32) -> Self {
        let n = ploidy * max_amplification;
        DiscreteAlleleFreqs {
            inner: (0..n + 1)
                .map(|m| AlleleFreq(m as f64 / n as f64))
                .collect_vec(),
        }
    }

    /// Return spectrum of possible allele frequencies given a ploidy.
    pub fn feasible(ploidy: u32) -> Self {
        Self::_feasible(ploidy, 1)
    }

    /// Return all frequencies except 0.0.
    pub fn not_absent(&self) -> Self {
        DiscreteAlleleFreqs {
            inner: self.inner[1..].to_owned(),
        }
    }

    pub fn absent() -> Self {
        DiscreteAlleleFreqs {
            inner: vec![AlleleFreq(0.0)],
        }
    }
}

impl Deref for DiscreteAlleleFreqs {
    type Target = Vec<AlleleFreq>;

    fn deref(&self) -> &Vec<AlleleFreq> {
        &self.inner
    }
}

/// An allele frequency range
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ContinuousAlleleFreqs {
    inner: Range<AlleleFreq>,
    pub left_exclusive: bool,
    pub right_exclusive: bool,
    /// offset to add when calculating the smallest observable value for a left-exclusive 0.0 bound
    zero_offset: NotNan<f64>,
}

impl ContinuousAlleleFreqs {
    pub fn absent() -> Self {
        Self::singleton(0.0)
    }

    pub fn singleton(value: f64) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(value)..AlleleFreq(value),
            left_exclusive: false,
            right_exclusive: false,
            zero_offset: NotNan::from(1.0),
        }
    }

    /// create a left- and right-inclusive allele frequency range
    pub fn inclusive(range: Range<f64>) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: false,
            right_exclusive: false,
            zero_offset: NotNan::from(1.0),
        }
    }

    /// create a left- and right-exclusive allele frequency range
    pub fn exclusive(range: Range<f64>) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: true,
            right_exclusive: true,
            zero_offset: NotNan::from(1.0),
        }
    }

    /// create a left-exclusive allele frequency range
    pub fn left_exclusive(range: Range<f64>) -> Self {
        if range.start == range.end {
            panic!("ContinuousAlleleFreqs::left_exclusive({}..{}) does not make sense with identical start and end point.", range.start, range.end);
        }
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: true,
            right_exclusive: false,
            zero_offset: NotNan::from(1.0),
        }
    }

    /// create a right-exclusive allele frequency range
    pub fn right_exclusive(range: Range<f64>) -> Self {
        if range.start == range.end {
            panic!("ContinuousAlleleFreqs::right_exclusive({}..{}) does not make sense with identical start and end point.", range.start, range.end);
        }
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: false,
            right_exclusive: true,
            zero_offset: NotNan::from(1.0),
        }
    }

    pub fn min_observations(mut self, min_observations: usize) -> Self {
        self.zero_offset = NotNan::from(min_observations as f64);

        self
    }

    pub fn is_singleton(&self) -> bool {
        self.start == self.end
    }

    pub fn observable_min(&self, n_obs: usize) -> AlleleFreq {
        if n_obs < 10 {
            self.start
        } else {
            let obs_count = Self::expected_observation_count(self.start, n_obs);
            let adjust_allelefreq = |obs_count: f64| AlleleFreq(obs_count.ceil() / n_obs as f64);

            if self.left_exclusive && obs_count % 1.0 == 0.0 {
                // We are left exclusive and need to find a supremum from the right.
                let offsets = if *self.start == 0.0 {
                    // The lower bound is zero, hence we apply first any given zero offset if
                    // possible.
                    vec![*self.zero_offset, 1.0, 0.0]
                } else {
                    vec![1.0, 0.0]
                };

                let adjusted_end = self.observable_max(n_obs);

                for offset in offsets {
                    let adjusted_obs_count = obs_count + offset;
                    let adjusted_start = adjust_allelefreq(adjusted_obs_count);
                    if *adjusted_start <= 1.0 && adjusted_start <= adjusted_end {
                        return adjusted_start;
                    }
                }
            }

            adjust_allelefreq(obs_count)
        }
    }

    pub fn observable_max(&self, n_obs: usize) -> AlleleFreq {
        assert!(
            *self.end != 0.0,
            "bug: observable_max may not be called if end=0.0."
        );
        if n_obs < 10 {
            self.end
        } else {
            let mut obs_count = Self::expected_observation_count(self.end, n_obs);
            if self.right_exclusive && obs_count % 1.0 == 0.0 {
                obs_count -= 1.0;
            }
            AlleleFreq(obs_count.floor() / n_obs as f64)
        }
    }

    fn expected_observation_count(freq: AlleleFreq, n_obs: usize) -> f64 {
        n_obs as f64 * *freq
    }
}

impl Default for ContinuousAlleleFreqs {
    fn default() -> Self {
        Self::absent()
    }
}

impl Ord for ContinuousAlleleFreqs {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.inner.start.cmp(&other.start) {
            Ordering::Equal => self.inner.end.cmp(&other.end),
            ord @ _ => ord,
        }
    }
}

impl PartialOrd for ContinuousAlleleFreqs {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Deref for ContinuousAlleleFreqs {
    type Target = Range<AlleleFreq>;

    fn deref(&self) -> &Range<AlleleFreq> {
        &self.inner
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumString, EnumIter, IntoStaticStr)]
pub enum VariantType {
    #[strum(serialize = "INS")]
    Insertion(Option<Range<u32>>),
    #[strum(serialize = "DEL")]
    Deletion(Option<Range<u32>>),
    #[strum(serialize = "SNV")]
    SNV,
    #[strum(serialize = "REF")]
    None, // site with no suggested alternative allele
}

impl From<&str> for VariantType {
    fn from(string: &str) -> VariantType {
        match string {
            "INS" => VariantType::Insertion(None),
            "DEL" => VariantType::Deletion(None),
            "SNV" => VariantType::SNV,
            "REF" => VariantType::None,
            _ => panic!("bug: given string does not describe a valid variant type"),
        }
    }
}

#[derive(Clone, Debug)]
pub enum Variant {
    Deletion(u32),
    Insertion(Vec<u8>),
    SNV(u8),
    None,
}

impl Variant {
    pub fn has_fragment_evidence(&self) -> bool {
        match self {
            &Variant::Deletion(_) => true,
            &Variant::Insertion(_) => true,
            &Variant::SNV(_) => false,
            &Variant::None => false,
        }
    }

    pub fn is_single_base(&self) -> bool {
        match self {
            &Variant::SNV(_) | &Variant::None => true,
            _ => false,
        }
    }

    pub fn is_snv(&self) -> bool {
        match self {
            &Variant::SNV(_) => true,
            _ => false,
        }
    }

    pub fn is_indel(&self) -> bool {
        match self {
            &Variant::Deletion(_) => true,
            &Variant::Insertion(_) => true,
            &Variant::SNV(_) => false,
            &Variant::None => false,
        }
    }

    pub fn is_type(&self, vartype: &VariantType) -> bool {
        match (self, vartype) {
            (&Variant::Deletion(l), &VariantType::Deletion(Some(ref range))) => {
                l >= range.start && l < range.end
            }
            (&Variant::Insertion(_), &VariantType::Insertion(Some(ref range))) => {
                self.len() >= range.start && self.len() < range.end
            }
            (&Variant::Deletion(_), &VariantType::Deletion(None)) => true,
            (&Variant::Insertion(_), &VariantType::Insertion(None)) => true,
            (&Variant::SNV(_), &VariantType::SNV) => true,
            (&Variant::None, &VariantType::None) => true,
            _ => false,
        }
    }

    pub fn end(&self, start: u32) -> u32 {
        match self {
            &Variant::Deletion(length) => start + length,
            &Variant::Insertion(_) => start + 1, // end of insertion is the next regular base
            &Variant::SNV(_) | &Variant::None => start,
        }
    }

    pub fn centerpoint(&self, start: u32) -> u32 {
        match self {
            &Variant::Deletion(length) => start + length / 2,
            &Variant::Insertion(_) => start, // end of insertion is the next regular base
            &Variant::SNV(_) | &Variant::None => start,
        }
    }

    pub fn len(&self) -> u32 {
        match self {
            &Variant::Deletion(l) => l,
            &Variant::Insertion(ref s) => s.len() as u32,
            &Variant::SNV(_) => 1,
            &Variant::None => 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::model::evidence::{Evidence, Observation};
    use crate::utils;

    use bio::stats::LogProb;

    pub fn observation(prob_mapping: LogProb, prob_alt: LogProb, prob_ref: LogProb) -> Observation {
        Observation::new(
            prob_mapping,
            prob_mapping.ln_one_minus_exp(),
            prob_alt,
            prob_ref,
            utils::max_prob(prob_ref, prob_alt),
            LogProb::ln_one(),
            LogProb::ln_one(),
            LogProb::ln_one(),
            true,
            true,
            Evidence::dummy_alignment(),
        )
    }

    // fn setup_pairwise_test<'a>(
    // ) -> PairCaller<ContinuousAlleleFreqs, DiscreteAlleleFreqs, priors::TumorNormalModel> {
    //     let insert_size = InsertSize {
    //         mean: 250.0,
    //         sd: 50.0,
    //     };
    //     let prior_model = priors::TumorNormalModel::new(2, 30.0, 1.0, 1.0, 3e9 as u64, Prob(0.001));
    //     let case_sample = Sample::new(
    //         bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
    //         true,
    //         AlignmentProperties::default(insert_size),
    //         LatentVariableModel::new(1.0),
    //         constants::PROB_ILLUMINA_INS,
    //         constants::PROB_ILLUMINA_DEL,
    //         Prob(0.0),
    //         Prob(0.0),
    //         10,
    //         500,
    //         &[],
    //     );
    //     let control_sample = Sample::new(
    //         bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
    //         true,
    //         AlignmentProperties::default(insert_size),
    //         LatentVariableModel::new(1.0),
    //         constants::PROB_ILLUMINA_INS,
    //         constants::PROB_ILLUMINA_DEL,
    //         Prob(0.0),
    //         Prob(0.0),
    //         10,
    //         500,
    //         &[],
    //     );
    //
    //     let model = PairCaller::new(case_sample, control_sample, prior_model);
    //
    //     model
    // }
    //
    //
    // /// scenario 1: same pileup -> germline call
    // #[test]
    // fn test_same_pileup() {
    //     let variant = Variant::Deletion(3);
    //     let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
    //     let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
    //     let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
    //     let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);
    //
    //     let mut observations = Vec::new();
    //     for _ in 0..5 {
    //         observations.push(observation(
    //             LogProb::ln_one(),
    //             LogProb::ln_one(),
    //             LogProb::ln_zero(),
    //         ));
    //     }
    //
    //     let model = setup_pairwise_test();
    //     let mut pileup = PairPileup::new(
    //         observations.clone(),
    //         observations.clone(),
    //         variant,
    //         &model.prior_model,
    //         model.case_sample.borrow().likelihood_model(),
    //         model.control_sample.borrow().likelihood_model(),
    //     );
    //
    //     // germline
    //     let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
    //     let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
    //     assert_relative_eq!(p_germline.exp(), 1.0);
    //     // somatic
    //     assert_relative_eq!(p_somatic.exp(), 0.0);
    //     assert_relative_eq!(
    //         p_germline.ln_add_exp(p_somatic).exp(),
    //         1.0,
    //         epsilon = 0.0000001
    //     );
    //     assert!(*pileup.map_allele_freqs().0 >= 0.97);
    // }
    //
    // /// scenario 2: empty control pileup -> somatic call
    // #[test]
    // fn test_empty_control_pileup() {
    //     let variant = Variant::Deletion(3);
    //     let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
    //     let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
    //     let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
    //     let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);
    //
    //     let mut observations = Vec::new();
    //     for _ in 0..5 {
    //         observations.push(observation(
    //             LogProb::ln_one(),
    //             LogProb::ln_one(),
    //             LogProb::ln_zero(),
    //         ));
    //     }
    //
    //     let model = setup_pairwise_test();
    //     let mut pileup = PairPileup::new(
    //         observations.clone(),
    //         vec![],
    //         variant,
    //         &model.prior_model,
    //         model.case_sample.borrow().likelihood_model(),
    //         model.control_sample.borrow().likelihood_model(),
    //     );
    //
    //     let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
    //     let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
    //     let (af_case, af_control) = pileup.map_allele_freqs();
    //     // we have no evidence for germline, but an allele frequency of 1 is most likely with a germline variant!
    //     assert!(p_germline > p_somatic);
    //     // germline < somatic
    //     assert_relative_eq!(
    //         p_germline.ln_add_exp(p_somatic).exp(),
    //         1.0,
    //         epsilon = 0.0000001
    //     );
    //     assert!(*af_case >= 0.97);
    //     assert!(*af_control >= 0.97);
    // }
    //
    // /// scenario 3: subclonal variant
    // #[test]
    // fn test_subclonal() {
    //     let variant = Variant::Deletion(3);
    //     let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
    //     let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
    //     let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
    //     let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);
    //
    //     let mut observations = Vec::new();
    //     for _ in 0..5 {
    //         observations.push(observation(
    //             LogProb::ln_one(),
    //             LogProb::ln_one(),
    //             LogProb::ln_zero(),
    //         ));
    //     }
    //     for _ in 0..50 {
    //         observations.push(observation(
    //             LogProb::ln_one(),
    //             LogProb::ln_zero(),
    //             LogProb::ln_one(),
    //         ));
    //     }
    //
    //     let model = setup_pairwise_test();
    //     let mut pileup = PairPileup::new(
    //         observations.clone(),
    //         vec![],
    //         variant,
    //         &model.prior_model,
    //         model.case_sample.borrow().likelihood_model(),
    //         model.control_sample.borrow().likelihood_model(),
    //     );
    //
    //     let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
    //     let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
    //     // somatic
    //     assert_relative_eq!(p_somatic.exp(), 0.9985, epsilon = 0.01);
    //     assert_relative_eq!(
    //         p_germline.ln_add_exp(p_somatic).exp(),
    //         1.0,
    //         epsilon = 0.0000001
    //     );
    //     assert_relative_eq!(*pileup.map_allele_freqs().0, 0.09, epsilon = 0.03);
    // }
    //
    // /// scenario 4: absent variant
    // #[test]
    // fn test_absent() {
    //     let variant = Variant::Deletion(3);
    //     let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
    //     let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
    //     let tumor_ref = ContinuousAlleleFreqs::inclusive(0.0..0.0);
    //     let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
    //     let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);
    //
    //     let mut observations = Vec::new();
    //     for _ in 0..10 {
    //         observations.push(observation(
    //             LogProb::ln_one(),
    //             LogProb::ln_zero(),
    //             LogProb::ln_one(),
    //         ));
    //     }
    //
    //     let model = setup_pairwise_test();
    //     let mut pileup = PairPileup::new(
    //         observations.clone(),
    //         observations.clone(),
    //         variant,
    //         &model.prior_model,
    //         model.case_sample.borrow().likelihood_model(),
    //         model.control_sample.borrow().likelihood_model(),
    //     );
    //
    //     let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
    //     let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
    //     let p_absent = pileup.posterior_prob(&tumor_ref, &normal_ref);
    //
    //     // germline
    //     assert_relative_eq!(p_germline.exp(), 0.0, epsilon = 0.02);
    //     // somatic
    //     assert_relative_eq!(p_somatic.exp(), 0.0, epsilon = 0.02);
    //     // absent
    //     assert_relative_eq!(p_absent.exp(), 1.0, epsilon = 0.02);
    // }
    //
    // #[test]
    // fn test_allele_freq_observable() {
    //     let afs = ContinuousAlleleFreqs::left_exclusive(0.0..0.5);
    //     assert_relative_eq!(*afs.observable_min(15), 0.0666666666666667);
    //     assert_relative_eq!(*afs.observable_max(15), 0.4666666666666667);
    //
    //     let afs = ContinuousAlleleFreqs::right_exclusive(0.0..0.5);
    //     assert_relative_eq!(*afs.observable_min(12), 0.0);
    //     assert_relative_eq!(*afs.observable_max(12), 0.4166666666666667);
    // }
}

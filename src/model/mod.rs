use std::cell::Cell;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Deref, Range};

use itertools::Itertools;
use ordered_float::NotNaN;

pub mod evidence;
pub mod likelihood;
pub mod priors;
pub mod sample;

use bio::stats::LogProb;

use model::evidence::Observation;
use model::sample::Sample;

pub type AlleleFreq = NotNaN<f64>;

#[allow(non_snake_case)]
pub fn AlleleFreq(af: f64) -> AlleleFreq {
    NotNaN::new(af).unwrap()
}

pub trait AlleleFreqs: Debug {}
impl AlleleFreqs for DiscreteAlleleFreqs {}
impl AlleleFreqs for ContinuousAlleleFreqs {}

#[derive(Debug, Clone)]
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
#[derive(Debug)]
pub struct ContinuousAlleleFreqs {
    inner: Range<AlleleFreq>,
    pub left_exclusive: bool,
    pub right_exclusive: bool,
}

impl ContinuousAlleleFreqs {
    /// create a left- and right-inclusive allele frequency range
    pub fn inclusive(range: Range<f64>) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: false,
            right_exclusive: false,
        }
    }

    /// create a left- and right-exclusive allele frequency range
    pub fn exclusive(range: Range<f64>) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: true,
            right_exclusive: true,
        }
    }

    /// create a left-exclusive allele frequency range
    pub fn left_exclusive(range: Range<f64>) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: true,
            right_exclusive: false,
        }
    }

    /// create a right-exclusive allele frequency range
    pub fn right_exclusive(range: Range<f64>) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: false,
            right_exclusive: true,
        }
    }
}

impl Deref for ContinuousAlleleFreqs {
    type Target = Range<AlleleFreq>;

    fn deref(&self) -> &Range<AlleleFreq> {
        &self.inner
    }
}

#[derive(Debug)]
pub enum VariantType {
    Insertion(Option<Range<u32>>),
    Deletion(Option<Range<u32>>),
    SNV,
    None, // site with no suggested alternative allele
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

/// Joint variant calling model, combining two latent variable models.
pub struct SingleCaller<A: AlleleFreqs, P: priors::Model<A>> {
    sample: RefCell<Sample>,
    prior_model: P,
    a: PhantomData<A>,
}

impl<A: AlleleFreqs, P: priors::Model<A>> SingleCaller<A, P> {
    /// Create new `JointModel`.
    ///
    /// # Arguments
    ///
    /// * `sample` - sample
    /// * `prior_model` - prior model
    pub fn new(sample: Sample, prior_model: P) -> Self {
        SingleCaller {
            sample: RefCell::new(sample),
            prior_model: prior_model,
            a: PhantomData,
        }
    }

    /// Calculate pileup and marginal probability for given variant.
    ///
    /// # Arguments
    ///
    /// * `chrom` - the chromosome of the variant
    /// * `start` - the starting position of the variant
    /// * `variant` - the variant
    ///
    /// # Returns
    /// The `SinglePileup`, or an error message.
    pub fn pileup(
        &self,
        chrom: &[u8],
        start: u32,
        variant: Variant,
        chrom_seq: &[u8],
    ) -> Result<SinglePileup<A, P>, Box<Error>> {
        let pileup = self
            .sample
            .borrow_mut()
            .extract_observations(start, &variant, chrom, chrom_seq)?;
        debug!("Obtained pileups ({} observations).", pileup.len());

        Ok(SinglePileup::new(
            pileup,
            variant,
            &self.prior_model,
            self.sample.borrow().likelihood_model(),
        ))
    }
}

/// Pileup of observations associated with marginal probability.
pub struct SinglePileup<'a, A, P>
where
    A: AlleleFreqs,
    P: 'a + priors::Model<A>,
{
    observations: Vec<Observation>,
    // we use Cell for marginal prob to be able to mutate the field without having mutable access to the whole pileup
    marginal_prob: Cell<Option<LogProb>>,
    prior_model: &'a P,
    sample_model: likelihood::LatentVariableModel,
    variant: Variant,
    a: PhantomData<A>,
}

impl<'a, A: AlleleFreqs, P: priors::Model<A>> SinglePileup<'a, A, P> {
    /// Create new pileup.
    fn new(
        observations: Vec<Observation>,
        variant: Variant,
        prior_model: &'a P,
        sample_model: likelihood::LatentVariableModel,
    ) -> Self {
        SinglePileup {
            observations: observations,
            marginal_prob: Cell::new(None),
            variant: variant,
            prior_model: prior_model,
            sample_model: sample_model,
            a: PhantomData,
        }
    }

    fn marginal_prob(&self) -> LogProb {
        if self.marginal_prob.get().is_none() {
            debug!("Calculating marginal probability.");

            let likelihood = |af: AlleleFreq| self.likelihood(af);
            let p =
                self.prior_model
                    .marginal_prob(&likelihood, &self.variant, self.observations.len());

            self.marginal_prob.set(Some(p));
        }

        self.marginal_prob.get().unwrap()
    }

    pub fn joint_prob(&self, af: &A) -> LogProb {
        let likelihood = |af: AlleleFreq| self.likelihood(af);

        let p =
            self.prior_model
                .joint_prob(af, &likelihood, &self.variant, self.observations.len());
        p
    }

    /// Calculate posterior probability of given allele frequencies.
    pub fn posterior_prob(&self, af: &A) -> LogProb {
        let p = self.joint_prob(af);
        let marginal = self.marginal_prob();
        let prob = p - marginal;
        prob
    }

    pub fn map_allele_freqs(&self) -> AlleleFreq {
        let likelihood = |af: AlleleFreq| self.likelihood(af);

        self.prior_model
            .map(&likelihood, &self.variant, self.observations.len())
    }

    fn likelihood(&self, af: AlleleFreq) -> LogProb {
        self.sample_model
            .likelihood_pileup(&self.observations, af, None)
    }

    pub fn observations(&self) -> &[Observation] {
        &self.observations
    }
}

/// Joint variant calling model, combining two latent variable models.
pub struct PairCaller<A: AlleleFreqs, B: AlleleFreqs, P: priors::PairModel<A, B>> {
    case_sample: RefCell<Sample>,
    control_sample: RefCell<Sample>,
    prior_model: P,
    a: PhantomData<A>,
    b: PhantomData<B>,
}

impl<A: AlleleFreqs, B: AlleleFreqs, P: priors::PairModel<A, B>> PairCaller<A, B, P> {
    /// Create new `PairCaller`.
    ///
    /// # Arguments
    ///
    /// * `case_sample` - case sample
    /// * `control_sample` - control sample
    /// * `prior_model` - prior model
    pub fn new(case_sample: Sample, control_sample: Sample, prior_model: P) -> Self {
        PairCaller {
            case_sample: RefCell::new(case_sample),
            control_sample: RefCell::new(control_sample),
            prior_model: prior_model,
            a: PhantomData,
            b: PhantomData,
        }
    }

    /// Calculate pileup and marginal probability for given variant.
    ///
    /// # Arguments
    ///
    /// * `chrom` - the chromosome of the variant
    /// * `start` - the starting position of the variant
    /// * `variant` - the variant
    ///
    /// # Returns
    /// The `PairPileup`, or an error message.
    pub fn pileup(
        &self,
        chrom: &[u8],
        start: u32,
        variant: Variant,
        chrom_seq: &[u8],
    ) -> Result<PairPileup<A, B, P>, Box<Error>> {
        debug!("Case pileup");
        let case_pileup = self
            .case_sample
            .borrow_mut()
            .extract_observations(start, &variant, chrom, chrom_seq)?;
        debug!("Control pileup");
        let control_pileup = self
            .control_sample
            .borrow_mut()
            .extract_observations(start, &variant, chrom, chrom_seq)?;

        debug!(
            "Obtained pileups (case: {} observations, control: {} observations).",
            case_pileup.len(),
            control_pileup.len()
        );
        Ok(PairPileup::new(
            case_pileup,
            control_pileup,
            variant,
            &self.prior_model,
            self.case_sample.borrow().likelihood_model(),
            self.control_sample.borrow().likelihood_model(),
        ))
    }
}

/// Pileup of observations associated with marginal probability.
pub struct PairPileup<'a, A, B, P: ?Sized>
where
    A: AlleleFreqs,
    B: AlleleFreqs,
    P: 'a + priors::PairModel<A, B>,
{
    case: Vec<Observation>,
    control: Vec<Observation>,
    // we use Cell for marginal prob to be able to mutate the field without having mutable access to the whole pileup
    marginal_prob: Cell<Option<LogProb>>,
    prior_model: &'a P,
    case_sample_model: likelihood::LatentVariableModel,
    control_sample_model: likelihood::LatentVariableModel,
    // these two caches are useful for PairModels where case and control are independent of each other and Events are just combinations of across the two axes, as e.g. in the single-cell-bulk model
    case_likelihood_cache: BTreeMap<AlleleFreq, LogProb>,
    control_likelihood_cache: BTreeMap<AlleleFreq, LogProb>,
    variant: Variant,
    a: PhantomData<A>,
    b: PhantomData<B>,
}

impl<'a, A: AlleleFreqs, B: AlleleFreqs, P: priors::PairModel<A, B>> PairPileup<'a, A, B, P> {
    /// Create new pileup.
    fn new(
        case: Vec<Observation>,
        control: Vec<Observation>,
        variant: Variant,
        prior_model: &'a P,
        case_sample_model: likelihood::LatentVariableModel,
        control_sample_model: likelihood::LatentVariableModel,
    ) -> Self {
        PairPileup {
            case: case,
            control: control,
            marginal_prob: Cell::new(None),
            variant: variant,
            prior_model: prior_model,
            case_sample_model: case_sample_model,
            control_sample_model: control_sample_model,
            case_likelihood_cache: BTreeMap::new(),
            control_likelihood_cache: BTreeMap::new(),
            a: PhantomData,
            b: PhantomData,
        }
    }

    fn marginal_prob(&mut self) -> LogProb {
        if self.marginal_prob.get().is_none() {
            debug!("Calculating marginal probability.");

            let p = self.prior_model.marginal_prob(self);
            debug!("Marginal probability: {}.", p.exp());

            self.marginal_prob.set(Some(p));
        }

        self.marginal_prob.get().unwrap()
    }

    #[cfg_attr(feature = "flame_it", flame)]
    pub fn joint_prob(&mut self, af_case: &A, af_control: &B) -> LogProb {
        let p = self.prior_model.joint_prob(af_case, af_control, self);
        p
    }

    /// Calculate posterior probability of given allele frequencies.
    pub fn posterior_prob(&mut self, af_case: &A, af_control: &B) -> LogProb {
        let p = self.joint_prob(af_case, af_control);
        let marginal = self.marginal_prob();
        let prob = p - marginal;
        prob
    }

    #[cfg_attr(feature = "flame_it", flame)]
    pub fn map_allele_freqs(&mut self) -> (AlleleFreq, AlleleFreq) {
        self.prior_model.map(self)
    }

    #[cfg_attr(feature = "flame_it", flame)]
    fn case_likelihood(&mut self, af_case: AlleleFreq, af_control: Option<AlleleFreq>) -> LogProb {
        // no af_control given, because case and control are independent
        if af_control.is_none() {
            // get likelihood if already cached
            if let Some(likelihood) = self.case_likelihood_cache.get(&af_case) {
                return *likelihood;
            }
            // compute otherwise
            let l = self
                .case_sample_model
                .likelihood_pileup(&self.case, af_case, None);
            // insert into cache before returning
            self.case_likelihood_cache.insert(af_case, l);
            l
        // cache cannot be used
        } else {
            self.case_sample_model
                .likelihood_pileup(&self.case, af_case, af_control)
        }
    }

    #[cfg_attr(feature = "flame_it", flame)]
    fn control_likelihood(
        &mut self,
        af_control: AlleleFreq,
        af_case: Option<AlleleFreq>,
    ) -> LogProb {
        // no af_control given, because case and control are independent
        if af_case.is_none() {
            // get likelihood if already cached
            if let Some(likelihood) = self.control_likelihood_cache.get(&af_control) {
                return *likelihood;
            }
            // compute otherwise
            let l = self
                .control_sample_model
                .likelihood_pileup(&self.control, af_control, None);
            // insert into cache before returning
            self.control_likelihood_cache.insert(af_control, l);
            l
        // cache cannot be used
        } else {
            self.control_sample_model
                .likelihood_pileup(&self.control, af_control, af_case)
        }
    }

    pub fn case_observations(&self) -> &[Observation] {
        &self.case
    }

    pub fn control_observations(&self) -> &[Observation] {
        &self.control
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use constants;
    use estimation::alignment_properties::{AlignmentProperties, InsertSize};
    use likelihood::LatentVariableModel;
    use model::evidence::{Evidence, Observation};
    use Sample;

    use bio::stats::{LogProb, Prob};
    use rust_htslib::bam;

    fn setup_pairwise_test<'a>(
) -> PairCaller<ContinuousAlleleFreqs, DiscreteAlleleFreqs, priors::TumorNormalModel> {
        let insert_size = InsertSize {
            mean: 250.0,
            sd: 50.0,
        };
        let prior_model = priors::TumorNormalModel::new(2, 30.0, 1.0, 1.0, 3e9 as u64, Prob(0.001));
        let case_sample = Sample::new(
            bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            true,
            AlignmentProperties::default(insert_size),
            LatentVariableModel::new(1.0),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            10,
        );
        let control_sample = Sample::new(
            bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            true,
            AlignmentProperties::default(insert_size),
            LatentVariableModel::new(1.0),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            10,
        );

        let model = PairCaller::new(case_sample, control_sample, prior_model);

        model
    }

    pub fn observation(prob_mapping: LogProb, prob_alt: LogProb, prob_ref: LogProb) -> Observation {
        Observation::new(
            prob_mapping,
            prob_alt,
            prob_ref,
            LogProb::ln_one(),
            Evidence::dummy_alignment(),
        )
    }

    /// scenario 1: same pileup -> germline call
    #[test]
    fn test_same_pileup() {
        let variant = Variant::Deletion(3);
        let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
        let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
        let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
        let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);

        let mut observations = Vec::new();
        for _ in 0..5 {
            observations.push(observation(
                LogProb::ln_one(),
                LogProb::ln_one(),
                LogProb::ln_zero(),
            ));
        }

        let model = setup_pairwise_test();
        let mut pileup = PairPileup::new(
            observations.clone(),
            observations.clone(),
            variant,
            &model.prior_model,
            model.case_sample.borrow().likelihood_model(),
            model.control_sample.borrow().likelihood_model(),
        );

        // germline
        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        assert_relative_eq!(p_germline.exp(), 1.0);
        // somatic
        assert_relative_eq!(p_somatic.exp(), 0.0);
        assert_relative_eq!(
            p_germline.ln_add_exp(p_somatic).exp(),
            1.0,
            epsilon = 0.0000001
        );
        assert!(*pileup.map_allele_freqs().0 >= 0.97);
    }

    /// scenario 2: empty control pileup -> somatic call
    #[test]
    fn test_empty_control_pileup() {
        let variant = Variant::Deletion(3);
        let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
        let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
        let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
        let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);

        let mut observations = Vec::new();
        for _ in 0..5 {
            observations.push(observation(
                LogProb::ln_one(),
                LogProb::ln_one(),
                LogProb::ln_zero(),
            ));
        }

        let model = setup_pairwise_test();
        let mut pileup = PairPileup::new(
            observations.clone(),
            vec![],
            variant,
            &model.prior_model,
            model.case_sample.borrow().likelihood_model(),
            model.control_sample.borrow().likelihood_model(),
        );

        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        let (af_case, af_control) = pileup.map_allele_freqs();
        // we have no evidence for germline, but an allele frequency of 1 is most likely with a germline variant!
        assert!(p_germline > p_somatic);
        // germline < somatic
        assert_relative_eq!(
            p_germline.ln_add_exp(p_somatic).exp(),
            1.0,
            epsilon = 0.0000001
        );
        assert!(*af_case >= 0.97);
        assert!(*af_control >= 0.97);
    }

    /// scenario 3: subclonal variant
    #[test]
    fn test_subclonal() {
        let variant = Variant::Deletion(3);
        let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
        let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
        let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
        let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);

        let mut observations = Vec::new();
        for _ in 0..5 {
            observations.push(observation(
                LogProb::ln_one(),
                LogProb::ln_one(),
                LogProb::ln_zero(),
            ));
        }
        for _ in 0..50 {
            observations.push(observation(
                LogProb::ln_one(),
                LogProb::ln_zero(),
                LogProb::ln_one(),
            ));
        }

        let model = setup_pairwise_test();
        let mut pileup = PairPileup::new(
            observations.clone(),
            vec![],
            variant,
            &model.prior_model,
            model.case_sample.borrow().likelihood_model(),
            model.control_sample.borrow().likelihood_model(),
        );

        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        // somatic
        assert_relative_eq!(p_somatic.exp(), 0.9985, epsilon = 0.01);
        assert_relative_eq!(
            p_germline.ln_add_exp(p_somatic).exp(),
            1.0,
            epsilon = 0.0000001
        );
        assert_relative_eq!(*pileup.map_allele_freqs().0, 0.09, epsilon = 0.03);
    }

    /// scenario 4: absent variant
    #[test]
    fn test_absent() {
        let variant = Variant::Deletion(3);
        let tumor_all = ContinuousAlleleFreqs::inclusive(0.0..1.0);
        let tumor_alt = ContinuousAlleleFreqs::left_exclusive(0.0..1.0);
        let tumor_ref = ContinuousAlleleFreqs::inclusive(0.0..0.0);
        let normal_alt = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5), AlleleFreq(1.0)]);
        let normal_ref = DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]);

        let mut observations = Vec::new();
        for _ in 0..10 {
            observations.push(observation(
                LogProb::ln_one(),
                LogProb::ln_zero(),
                LogProb::ln_one(),
            ));
        }

        let model = setup_pairwise_test();
        let mut pileup = PairPileup::new(
            observations.clone(),
            observations.clone(),
            variant,
            &model.prior_model,
            model.case_sample.borrow().likelihood_model(),
            model.control_sample.borrow().likelihood_model(),
        );

        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        let p_absent = pileup.posterior_prob(&tumor_ref, &normal_ref);

        // germline
        assert_relative_eq!(p_germline.exp(), 0.0, epsilon = 0.02);
        // somatic
        assert_relative_eq!(p_somatic.exp(), 0.0, epsilon = 0.02);
        // absent
        assert_relative_eq!(p_absent.exp(), 1.0, epsilon = 0.02);
    }
}

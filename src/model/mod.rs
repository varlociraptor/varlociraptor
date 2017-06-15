use std::marker::PhantomData;
use std::cell::Cell;
use std::ops::{Range, Deref};
use std::fmt::Debug;
use std::error::Error;
use std::cell::RefCell;

use ordered_float::NotNaN;

pub mod likelihood;
pub mod priors;
pub mod sample;
pub mod evidence;

use bio::stats::LogProb;

use model::sample::{Sample, Observation};


pub type AlleleFreq = NotNaN<f64>;
pub type DiscreteAlleleFreqs = Vec<AlleleFreq>;


#[allow(non_snake_case)]
pub fn AlleleFreq(af: f64) -> AlleleFreq {
    NotNaN::new(af).unwrap()
}


pub trait AlleleFreqs: Debug {}
impl AlleleFreqs for DiscreteAlleleFreqs {}
impl AlleleFreqs for ContinuousAlleleFreqs {}

/// An allele frequency range
#[derive(Debug)]
pub struct ContinuousAlleleFreqs {
    inner: Range<AlleleFreq>,
    pub left_exclusive: bool,
    pub right_exclusive: bool
}

impl ContinuousAlleleFreqs {
    /// create a left- and right-inclusive allele frequency range
    pub fn inclusive(range: Range<f64>) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: false,
            right_exclusive: false
        }
    }

    /// create a left- and right-exclusive allele frequency range
    pub fn exclusive(range: Range<f64>) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: true,
            right_exclusive: true
        }
    }

    /// create a left-exclusive allele frequency range
    pub fn left_exclusive(range: Range<f64>) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: true,
            right_exclusive: false
        }
    }

    /// create a right-exclusive allele frequency range
    pub fn right_exclusive(range: Range<f64>) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(range.start)..AlleleFreq(range.end),
            left_exclusive: false,
            right_exclusive: true
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
    SNV
}


#[derive(Clone, Debug)]
pub enum Variant {
    Deletion(u32),
    Insertion(Vec<u8>),
    SNV(u8)
}


impl Variant {
    pub fn has_fragment_evidence(&self) -> bool {
        match self {
            &Variant::Deletion(_)  => true,
            &Variant::Insertion(_) => true,
            &Variant::SNV(_)       => false
        }
    }

    pub fn is_indel(&self) -> bool {
        match self {
            &Variant::Deletion(_)  => true,
            &Variant::Insertion(_) => true,
            &Variant::SNV(_)       => false
        }
    }

    pub fn is_type(&self, vartype: &VariantType) -> bool {
        match (self, vartype) {
            (&Variant::Deletion(l), &VariantType::Deletion(Some(ref range))) => {
                l >= range.start && l < range.end
            },
            (&Variant::Insertion(_), &VariantType::Insertion(Some(ref range))) => {
                self.len() >= range.start && self.len() < range.end
            },
            (&Variant::Deletion(_), &VariantType::Deletion(None)) => true,
            (&Variant::Insertion(_), &VariantType::Insertion(None)) => true,
            (&Variant::SNV(_), &VariantType::SNV) => true,
            _ => false
        }
    }

    pub fn len(&self) -> u32 {
        match self {
            &Variant::Deletion(l)      => l,
            &Variant::Insertion(ref s) => s.len() as u32,
            &Variant::SNV(_)           => 1
        }
    }
}


/// Joint variant calling model, combining two latent variable models.
pub struct SingleCaller<A: AlleleFreqs, P: priors::Model<A>> {
    sample: RefCell<Sample>,
    prior_model: P,
    a: PhantomData<A>
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
            a: PhantomData
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
    pub fn pileup(&self, chrom: &[u8], start: u32, variant: Variant, chrom_seq: &[u8]) -> Result<SinglePileup<A, P>, Box<Error>> {
        let pileup = try!(self.sample.borrow_mut().extract_observations(chrom, start, &variant, chrom_seq));
        debug!("Obtained pileups ({} observations).", pileup.len());
        Ok(SinglePileup::new(
            pileup,
            variant,
            &self.prior_model,
            self.sample.borrow().likelihood_model()
        ))
    }
}


/// Pileup of observations associated with marginal probability.
pub struct SinglePileup<'a, A, P> where
    A: AlleleFreqs,
    P: 'a + priors::Model<A>
{
    observations: Vec<Observation>,
    // we use Cell for marginal prob to be able to mutate the field without having mutable access to the whole pileup
    marginal_prob: Cell<Option<LogProb>>,
    prior_model: &'a P,
    sample_model: likelihood::LatentVariableModel,
    variant: Variant,
    a: PhantomData<A>
}


impl<'a, A: AlleleFreqs, P: priors::Model<A>> SinglePileup<'a, A, P> {
    /// Create new pileup.
    fn new(
        observations: Vec<Observation>,
        variant: Variant,
        prior_model: &'a P,
        sample_model: likelihood::LatentVariableModel
    ) -> Self {
        SinglePileup {
            observations: observations,
            marginal_prob: Cell::new(None),
            variant: variant,
            prior_model: prior_model,
            sample_model: sample_model,
            a: PhantomData
        }
    }

    fn marginal_prob(&self) -> LogProb {
        if self.marginal_prob.get().is_none() {
            debug!("Calculating marginal probability.");

            let likelihood = |af: AlleleFreq| {
                self.likelihood(af)
            };
            let p = self.prior_model.marginal_prob(&likelihood, &self.variant, self.observations.len());

            self.marginal_prob.set(Some(p));
        }

        self.marginal_prob.get().unwrap()
    }

    fn joint_prob(&self, af: &A) -> LogProb {
        let likelihood = |af: AlleleFreq| {
            self.likelihood(af)
        };

        let p = self.prior_model.joint_prob(af, &likelihood, &self.variant, self.observations.len());
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
        let likelihood = |af: AlleleFreq| {
            self.likelihood(af)
        };

        self.prior_model.map(&likelihood, &self.variant, self.observations.len())
    }

    fn likelihood(&self, af: AlleleFreq) -> LogProb {
        self.sample_model.likelihood_pileup(&self.observations, *af, None)
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
    b: PhantomData<B>
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
            b: PhantomData
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
    pub fn pileup(&self, chrom: &[u8], start: u32, variant: Variant, chrom_seq: &[u8]) -> Result<PairPileup<A, B, P>, Box<Error>> {
        debug!("Case pileup");
        let case_pileup = try!(self.case_sample.borrow_mut().extract_observations(chrom, start, &variant, chrom_seq));
        debug!("Control pileup");
        let control_pileup = try!(self.control_sample.borrow_mut().extract_observations(chrom, start, &variant, chrom_seq));
        debug!("Obtained pileups (case: {} observations, control: {} observations).", case_pileup.len(), control_pileup.len());
        Ok(PairPileup::new(
            case_pileup,
            control_pileup,
            variant,
            &self.prior_model,
            self.case_sample.borrow().likelihood_model(),
            self.control_sample.borrow().likelihood_model()
        ))
    }
}


/// Pileup of observations associated with marginal probability.
pub struct PairPileup<'a, A, B, P> where
    A: AlleleFreqs,
    B: AlleleFreqs,
    P: 'a + priors::PairModel<A, B>
{
    case: Vec<Observation>,
    control: Vec<Observation>,
    // we use Cell for marginal prob to be able to mutate the field without having mutable access to the whole pileup
    marginal_prob: Cell<Option<LogProb>>,
    prior_model: &'a P,
    case_sample_model: likelihood::LatentVariableModel,
    control_sample_model: likelihood::LatentVariableModel,
    variant: Variant,
    a: PhantomData<A>,
    b: PhantomData<B>
}


impl<'a, A: AlleleFreqs, B: AlleleFreqs, P: priors::PairModel<A, B>> PairPileup<'a, A, B, P> {
    /// Create new pileup.
    fn new(
        case: Vec<Observation>,
        control: Vec<Observation>,
        variant: Variant,
        prior_model: &'a P,
        case_sample_model: likelihood::LatentVariableModel,
        control_sample_model: likelihood::LatentVariableModel
    ) -> Self {
        PairPileup {
            case: case,
            control: control,
            marginal_prob: Cell::new(None),
            variant: variant,
            prior_model: prior_model,
            case_sample_model: case_sample_model,
            control_sample_model: control_sample_model,
            a: PhantomData,
            b: PhantomData
        }
    }

    fn marginal_prob(&self) -> LogProb {
        if self.marginal_prob.get().is_none() {
            debug!("Calculating marginal probability.");

            let case_likelihood = |af_case: AlleleFreq, af_control: Option<AlleleFreq>| {
                self.case_likelihood(af_case, af_control)
            };
            let control_likelihood = |af_control: AlleleFreq, af_case: Option<AlleleFreq>| {
                self.control_likelihood(af_control, af_case)
            };
            let p = self.prior_model.marginal_prob(&case_likelihood, &control_likelihood, self.variant, self.case.len(), self.control.len());
            debug!("Marginal probability: {}.", p.exp());

            self.marginal_prob.set(Some(p));
        }

        self.marginal_prob.get().unwrap()
    }

    fn joint_prob(&self, af_case: &A, af_control: &B) -> LogProb {
        let case_likelihood = |af_case: AlleleFreq, af_control: Option<AlleleFreq>| {
            self.case_likelihood(af_case, af_control)
        };
        let control_likelihood = |af_control: AlleleFreq, af_case: Option<AlleleFreq>| {
            self.control_likelihood(af_control, af_case)
        };

        let p = self.prior_model.joint_prob(af_case, af_control, &case_likelihood, &control_likelihood, &self.variant, self.case.len(), self.control.len());
        debug!("Pr(D, f_case={:?}, f_control={:?}) = {}", af_case, af_control, p.exp());
        p
    }

    /// Calculate posterior probability of given allele frequencies.
    pub fn posterior_prob(&self, af_case: &A, af_control: &B) -> LogProb {
        debug!("Calculating posterior probability. for case={:?} and control={:?}", af_case, af_control);
        let p = self.joint_prob(af_case, af_control);
        let marginal = self.marginal_prob();
        let prob = p - marginal;
        debug!("Posterior probability: {}.", prob.exp());
        prob
    }

    pub fn map_allele_freqs(&self) -> (AlleleFreq, AlleleFreq) {
        let case_likelihood = |af_case: AlleleFreq, af_control: Option<AlleleFreq>| {
            self.case_likelihood(af_case, af_control)
        };
        let control_likelihood = |af_control: AlleleFreq, af_case: Option<AlleleFreq>| {
            self.control_likelihood(af_control, af_case)
        };

        self.prior_model.map(&case_likelihood, &control_likelihood, &self.variant, self.case.len(), self.control.len())
    }

    fn case_likelihood(&self, af_case: AlleleFreq, af_control: Option<AlleleFreq>) -> LogProb {
        self.case_sample_model.likelihood_pileup(&self.case, *af_case, af_control.map(|af| *af))
    }

    fn control_likelihood(&self, af_control: AlleleFreq, af_case: Option<AlleleFreq>) -> LogProb {
        self.control_sample_model.likelihood_pileup(&self.control, *af_control, af_case.map(|af| *af))
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
    use likelihood::LatentVariableModel;
    use Sample;
    use InsertSize;
    use model::sample::{Observation, Evidence};
    use constants;

    use rust_htslib::bam;
    use bio::stats::{LogProb, Prob};
    #[cfg(feature="flame_it")]
    use std::fs::File;
    #[cfg(feature="flame_it")]
    use flame;
    use csv;
    use itertools::Itertools;

    fn setup_pairwise_test<'a>() -> PairCaller<ContinuousAlleleFreqs, DiscreteAlleleFreqs, priors::TumorNormalModel> {
        let insert_size = InsertSize{ mean: 250.0, sd: 50.0 };
        let prior_model = priors::TumorNormalModel::new(2, 30.0, 1.0, 1.0, 3e9 as u64, Prob(0.001));
        let case_sample = Sample::new(
            bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            true,
            false,
            insert_size,
            LatentVariableModel::new(1.0),
            Prob(0.0),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            20,
            10
        );
        let control_sample = Sample::new(
            bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            true,
            false,
            insert_size,
            LatentVariableModel::new(1.0),
            Prob(0.0),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            20,
            10
        );

        let model = PairCaller::new(
            case_sample,
            control_sample,
            prior_model
        );

        model
    }

    /// scenario 1: same pileup -> germline call
    #[test]
    fn test_same_pileup() {
        let variant = Variant::Deletion(3);
        let tumor_all = ContinuousAlleleFreqs::inclusive( 0.0..1.0 );
        let tumor_alt = ContinuousAlleleFreqs::inclusive( 0.0..1.0 ); // TODO: should this be left_exclusive() instead?
        let normal_alt = vec![AlleleFreq(0.5), AlleleFreq(1.0)];
        let normal_ref = vec![AlleleFreq(0.0)];

        let mut observations = Vec::new();
        for _ in 0..5 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_one(),
                prob_ref: LogProb::ln_zero(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::dummy_alignment()
            });
        }

        let model = setup_pairwise_test();
        let pileup = PairPileup::new(
            observations.clone(), observations.clone(),
            variant,
            &model.prior_model,
            model.case_sample.borrow().likelihood_model(),
            model.control_sample.borrow().likelihood_model()
        );

        // germline
        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        assert_relative_eq!(p_germline.exp(), 1.0);
        // somatic
        assert_relative_eq!(p_somatic.exp(), 0.0);
        assert_relative_eq!(p_germline.ln_add_exp(p_somatic).exp(), 1.0, epsilon=0.0000001);
        assert!(*pileup.map_allele_freqs().0 >= 0.97);
    }

    /// scenario 2: empty control pileup -> somatic call
    #[test]
    fn test_empty_control_pileup() {
        let variant = Variant::Deletion(3);
        let tumor_all = ContinuousAlleleFreqs::inclusive( 0.0..1.0 );
        let tumor_alt = ContinuousAlleleFreqs::left_exclusive( 0.0..1.0 );
        let normal_alt = vec![AlleleFreq(0.5), AlleleFreq(1.0)];
        let normal_ref = vec![AlleleFreq(0.0)];

        let mut observations = Vec::new();
        for _ in 0..5 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_one(),
                prob_ref: LogProb::ln_zero(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::dummy_alignment()
            });
        }

        let model = setup_pairwise_test();
        let pileup = PairPileup::new(
            observations.clone(), vec![],
            variant,
            &model.prior_model,
            model.case_sample.borrow().likelihood_model(),
            model.control_sample.borrow().likelihood_model()
        );

        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        let (af_case, af_control) = pileup.map_allele_freqs();
        // we have no evidence for germline, but an allele frequency of 1 is most likely with a germline variant!
        assert!(p_germline > p_somatic);
        // germline < somatic
        assert_relative_eq!(p_germline.ln_add_exp(p_somatic).exp(), 1.0, epsilon=0.0000001);
        assert!(*af_case >= 0.97);
        assert!(*af_control >= 0.97);
    }

    /// scenario 3: subclonal variant
    #[test]
    fn test_subclonal() {
        let variant = Variant::Deletion(3);
        let tumor_all = ContinuousAlleleFreqs::inclusive( 0.0..1.0 );
        let tumor_alt = ContinuousAlleleFreqs::left_exclusive( 0.0..1.0 );
        let normal_alt = vec![AlleleFreq(0.5), AlleleFreq(1.0)];
        let normal_ref = vec![AlleleFreq(0.0)];

        let mut observations = Vec::new();
        for _ in 0..5 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_one(),
                prob_ref: LogProb::ln_zero(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::dummy_alignment()
            });
        }
        for _ in 0..50 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_zero(),
                prob_ref: LogProb::ln_one(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::dummy_alignment()
            });
        }

        let model = setup_pairwise_test();
        let pileup = PairPileup::new(
            observations.clone(), vec![],
            variant,
            &model.prior_model,
            model.case_sample.borrow().likelihood_model(),
            model.control_sample.borrow().likelihood_model()
        );

        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        // somatic
        assert_relative_eq!(p_somatic.exp(), 0.9985, epsilon=0.01);
        assert_relative_eq!(p_germline.ln_add_exp(p_somatic).exp(), 1.0, epsilon=0.0000001);
        assert_relative_eq!(*pileup.map_allele_freqs().0, 0.09, epsilon=0.03);
    }

    /// scenario 4: absent variant
    #[test]
    fn test_absent() {
        let variant = Variant::Deletion(3);
        let tumor_all = ContinuousAlleleFreqs::inclusive( 0.0..1.0 );
        let tumor_alt = ContinuousAlleleFreqs::left_exclusive( 0.0..1.0 );
        let tumor_ref = ContinuousAlleleFreqs::inclusive( 0.0..0.0 );
        let normal_alt = vec![AlleleFreq(0.5), AlleleFreq(1.0)];
        let normal_ref = vec![AlleleFreq(0.0)];

        let mut observations = Vec::new();
        for _ in 0..10 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_zero(),
                prob_ref: LogProb::ln_one(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::dummy_alignment()
            });
        }

        let model = setup_pairwise_test();
        let pileup = PairPileup::new(
            observations.clone(), observations.clone(),
            variant,
            &model.prior_model,
            model.case_sample.borrow().likelihood_model(),
            model.control_sample.borrow().likelihood_model()
        );

        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        let p_absent = pileup.posterior_prob(&tumor_ref, &normal_ref);

        // germline
        assert_relative_eq!(p_germline.exp(), 0.0, epsilon=0.02);
        // somatic
        assert_relative_eq!(p_somatic.exp(), 0.0, epsilon=0.02);
        // absent
        assert_relative_eq!(p_absent.exp(), 1.0, epsilon=0.02);
    }

    #[test]
    fn test_example1() {
        let variant = Variant::Insertion(b"AC".to_vec());
        let case_obs = vec![Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.507675873696745), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.0031672882261573254), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.507675873696745), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.0025150465111820103), prob_alt: LogProb::ln_one(), prob_ref: LogProb::ln_zero(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.0031672882261573254), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.507675873696745), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.0005013128699288086), prob_alt: LogProb::ln_one(), prob_ref: LogProb::ln_zero(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.007974998278512672), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.00000000010000000000499996), prob_alt: LogProb::ln_one(), prob_ref: LogProb::ln_zero(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.0031723009285603327), prob_alt: LogProb(-111.18254428986242), prob_ref: LogProb(-109.23587762319576), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.003994511005101995), prob_alt: LogProb(-113.14698873430689), prob_ref: LogProb(-111.18254428986242), prob_mismapped: LogProb::ln_one() }];
        let control_obs = vec![Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.0010005003335835337), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.0010005003335835337), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.00000000010000000000499996), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(), prob_mapping: LogProb(-0.00025122019630215495), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(),prob_mapping: LogProb(-0.0005013128699288086), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(),prob_mapping: LogProb(-0.0001585018800054507), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(),prob_mapping: LogProb(-0.0006311564818346603), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(),prob_mapping: LogProb(-0.00000000010000000000499996), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(),prob_mapping: LogProb(-0.0005013128699288086), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }, Observation { evidence: Evidence::dummy_alignment(),prob_mapping: LogProb(-0.003989017266406586), prob_alt: LogProb::ln_zero(), prob_ref: LogProb::ln_one(), prob_mismapped: LogProb::ln_one() }];

        let insert_size = InsertSize{ mean: 112.0, sd: 15.0 };
        let prior_model = priors::TumorNormalModel::new(2, 40000.0, 0.5, 0.5, 3e9 as u64, Prob(1.25E-4));
        let case_sample = Sample::new(
            bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            true,
            false,
            insert_size,
            LatentVariableModel::new(0.75),
            Prob(0.0),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            20,
            10
        );
        let control_sample = Sample::new(
            bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            true,
            false,
            insert_size,
            LatentVariableModel::new(1.0),
            Prob(0.0),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            20,
            10
        );

        let model = PairCaller::new(
            case_sample,
            control_sample,
            prior_model
        );

        let tumor_all = ContinuousAlleleFreqs::inclusive( 0.0..1.0 );
        let tumor_alt = ContinuousAlleleFreqs::inclusive( 0.05..1.0 ); // TODO: should this be left_exclusive() instead?
        let tumor_ref = ContinuousAlleleFreqs::inclusive( 0.0..0.001 ); // TODO: should this be (0.0)..(0.0) instead?
        let normal_alt = vec![AlleleFreq(0.5), AlleleFreq(1.0)];
        let normal_ref = vec![AlleleFreq(0.0)];

        let pileup = PairPileup::new(case_obs, control_obs, variant, &model.prior_model, model.case_sample.borrow().likelihood_model(), model.control_sample.borrow().likelihood_model());

        let p_absent = pileup.posterior_prob(&tumor_ref, &normal_ref);
        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        println!("{} {} {}", p_absent.exp(), p_somatic.exp(), p_germline.exp());

        assert!(p_somatic.exp() > 0.9);
        assert!(p_somatic > p_germline);
        assert!(p_somatic > p_absent);
    }

    fn setup_example(path: &str, deletion_factor: f64, insertion_factor: f64) -> (Vec<Observation>, Vec<Observation>, PairCaller<ContinuousAlleleFreqs, DiscreteAlleleFreqs, priors::TumorNormalModel>) {
        let mut reader = csv::Reader::from_file(path).expect("error reading example").delimiter(b'\t');
        let obs = reader.decode().collect::<Result<Vec<(String, u32, u32, String, Observation)>, _>>().unwrap();
        let groups = obs.into_iter().group_by(|&(_, _, _, ref sample, _)| {
            sample == "case"
        });
        let mut group_iter = groups.into_iter();
        let case_obs = group_iter.next().unwrap().1.into_iter().map(|(_, _, _, _, obs)| obs).collect_vec();
        let control_obs = if let Some(o) = group_iter.next() {
            o.1.into_iter().map(|(_, _, _, _, obs)| obs).collect_vec()
        } else {
            vec![]
        };

        let insert_size = InsertSize{ mean: 312.0, sd: 15.0 };
        let case_sample = Sample::new(
            bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            true,
            false,
            insert_size,
            LatentVariableModel::new(0.75),
            Prob(0.00001),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            20,
            10
        );
        let control_sample = Sample::new(
            bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            true,
            false,
            insert_size,
            LatentVariableModel::new(1.0),
            Prob(0.00001),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            20,
            10
        );


        let prior_model = priors::TumorNormalModel::new(2, 40000.0, deletion_factor, insertion_factor, 3e9 as u64, Prob(1.25E-4));
        let model = PairCaller::new(
            case_sample,
            control_sample,
            prior_model
        );

        (case_obs, control_obs, model)
    }

    #[allow(dead_code)]
    fn setup_example_flat(path: &str) -> (Vec<Observation>, Vec<Observation>, PairCaller<ContinuousAlleleFreqs, DiscreteAlleleFreqs, priors::FlatTumorNormalModel>) {
        let mut reader = csv::Reader::from_file(path).expect("error reading example").delimiter(b'\t');
        let obs = reader.decode().collect::<Result<Vec<(String, u32, u32, String, Observation)>, _>>().unwrap();
        let groups = obs.into_iter().group_by(|&(_, _, _, ref sample, _)| {
            sample == "case"
        });
        let mut group_iter = groups.into_iter();
        let case_obs = group_iter.next().unwrap().1.into_iter().map(|(_, _, _, _, obs)| obs).collect_vec();
        let control_obs = if let Some(o) = group_iter.next() {
            o.1.into_iter().map(|(_, _, _, _, obs)| obs).collect_vec()
        } else {
            vec![]
        };

        let insert_size = InsertSize{ mean: 312.0, sd: 15.0 };
        let case_sample = Sample::new(
            bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            true,
            false,
            insert_size,
            LatentVariableModel::new(0.75),
            Prob(0.00001),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            20,
            10
        );
        let control_sample = Sample::new(
            bam::IndexedReader::from_path(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            true,
            false,
            insert_size,
            LatentVariableModel::new(1.0),
            Prob(0.00001),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            20,
            10
        );


        let prior_model = priors::FlatTumorNormalModel::new(2);
        let model = PairCaller::new(
            case_sample,
            control_sample,
            prior_model
        );

        (case_obs, control_obs, model)
    }

    /// Test example where variant is subclonal and mapping quality is good
    #[test]
    fn test_example2() {
        let (case_obs, control_obs, model) = setup_example("tests/example2.obs.txt", 0.01, 0.03);

        let tumor_all = ContinuousAlleleFreqs::inclusive( 0.0..1.0 );
        let tumor_alt = ContinuousAlleleFreqs::inclusive( 0.05..1.0 ); // TODO: should this be left_exclusive() instead?
        let normal_alt = vec![AlleleFreq(0.5), AlleleFreq(1.0)];
        let normal_ref = vec![AlleleFreq(0.0)];

        let variant = Variant::Insertion(b"CTAT".to_vec());
        let pileup = PairPileup::new(case_obs, control_obs, variant, &model.prior_model, model.case_sample.borrow().likelihood_model(), model.control_sample.borrow().likelihood_model());

        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        let (af_case, _) = pileup.map_allele_freqs();
        println!("{} {}", p_somatic.exp(), p_germline.exp());

        assert!(p_somatic > p_germline);
        assert_relative_eq!(*af_case, 0.1, epsilon = 0.01);
    }

    /// Test example where variant is subclonal, but mapping quality is too bad
    #[test]
    fn test_example3() {
        let (case_obs, control_obs, model) = setup_example("tests/example3.obs.txt", 0.01, 0.03);

        let tumor_all = ContinuousAlleleFreqs::inclusive( 0.0..1.0 );
        let tumor_alt = ContinuousAlleleFreqs::inclusive( 0.05..1.0 ); // TODO: should this be left_exclusive() instead?
        let normal_alt = vec![AlleleFreq(0.5), AlleleFreq(1.0)];
        let normal_ref = vec![AlleleFreq(0.0)];

        let variant = Variant::Insertion(b"C".to_vec());
        let pileup = PairPileup::new(case_obs, control_obs, variant, &model.prior_model, model.case_sample.borrow().likelihood_model(), model.control_sample.borrow().likelihood_model());

        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        let (af_case, af_control) = pileup.map_allele_freqs();
        println!("{} {} {} {}", p_somatic.exp(), p_germline.exp(), af_case, af_control);
        assert!(p_somatic >= p_germline);
        assert!(*af_case <= 0.1);
        assert_relative_eq!(*af_control, 0.0);
    }

    /// Test example where variant is subclonal, but mapping quality is too bad
    #[test]
    fn test_example4() {
        let (case_obs, control_obs, model) = setup_example("tests/example4.obs.txt", 0.5, 0.5);
        /*for obs in case_obs.iter_mut() {
            obs.prob_mapping = LogProb(0.999f64.ln());
        }*/

        let tumor_all = ContinuousAlleleFreqs::inclusive( 0.0..1.0 );
        let tumor_alt = ContinuousAlleleFreqs::left_exclusive( 0.0..1.0 );
        let tumor_ref = ContinuousAlleleFreqs::inclusive( 0.0..0.0 );
        let normal_alt = vec![AlleleFreq(0.5), AlleleFreq(1.0)];
        let normal_ref = vec![AlleleFreq(0.0)];
        //let lh = model.case_sample().likelihood_model().likelihood_pileup(&case_obs[..100], 0.07, 0.0);
        let variant = Variant::Insertion(b"C".to_vec());
        /*let p_absent = model.joint_prob(&case_obs, &control_obs, &tumor_ref, &normal_ref, variant);
        let p_somatic = model.joint_prob(&case_obs, &control_obs, &tumor_alt, &normal_ref, variant);
        let p_germline = model.joint_prob(&case_obs, &control_obs, &tumor_all, &normal_alt, variant);
        let p_marginal = model.marginal_prob(&case_obs, &control_obs, variant);
        let norm = LogProb::ln_sum_exp(&[p_absent, p_somatic, p_germline]);

        println!("marginal={:e}, absent={}, somatic={}, germline={:e}, sum={:e}", p_marginal.exp(), (p_absent - norm).exp(), (p_somatic - norm).exp(), (p_germline - norm).exp(), norm.exp());
        assert!(false);*/

        let pileup = PairPileup::new(case_obs.to_owned(), control_obs, variant, &model.prior_model, model.case_sample.borrow().likelihood_model(), model.control_sample.borrow().likelihood_model());

        let p_somatic = pileup.posterior_prob(&tumor_alt, &normal_ref);
        let p_germline = pileup.posterior_prob(&tumor_all, &normal_alt);
        let p_absent = pileup.posterior_prob(&tumor_ref, &normal_ref);
        let (af_case, af_control) = pileup.map_allele_freqs();
        println!("{} {} {} {} {}", p_somatic.exp(), p_germline.exp(), p_absent.exp(), af_case, af_control);
        assert!(p_somatic >= p_germline);
        assert!(*af_case <= 0.1);
        assert_relative_eq!(*af_control, 0.0);
    }
}

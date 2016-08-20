use std::marker::PhantomData;
use std::cell::Cell;
use std::ops::Range;
use std::fmt::Debug;

use ordered_float::NotNaN;

pub mod likelihood;
pub mod priors;
pub mod sample;

use bio::stats::LogProb;

use model::sample::{Sample, Observation};


pub type AlleleFreq = NotNaN<f64>;
pub type DiscreteAlleleFreqs = Vec<AlleleFreq>;
pub type ContinousAlleleFreqs = Range<AlleleFreq>;


#[allow(non_snake_case)]
pub fn AlleleFreq(af: f64) -> AlleleFreq {
    NotNaN::new(af).unwrap()
}


pub trait AlleleFreqs: Debug {}
impl AlleleFreqs for DiscreteAlleleFreqs {}
impl AlleleFreqs for ContinousAlleleFreqs {}


#[derive(Copy, Clone)]
pub enum Variant {
    Deletion(u32),
    Insertion(u32)
}


pub trait JointModel<A: AlleleFreqs, B: AlleleFreqs, P: priors::PairModel<A, B>> {
    /// Return case sample.
    fn case_sample(&self) -> &Sample;

    /// Return control sample.
    fn control_sample(&self) -> &Sample;

    /// Return case sample.
    fn case_sample_mut(&mut self) -> &mut Sample;

    /// Return control sample.
    fn control_sample_mut(&mut self) -> &mut Sample;

    fn prior_model(&self) -> &P;

    /// Calculate joint probability of given pileups and allele frequencies.
    fn joint_prob(&self, case_pileup: &[Observation], control_pileup: &[Observation], af_case: &A, af_control: &B, variant: Variant) -> LogProb;

    /// Calculate marginal probability of given pileups.
    fn marginal_prob(&self, case_pileup: &[Observation], control_pileup: &[Observation], variant: Variant) -> LogProb {
        let af_case = self.prior_model().allele_freqs_case();
        let af_control = self.prior_model().allele_freqs_control();
        self.joint_prob(case_pileup, control_pileup, af_case, af_control, variant)
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
    /// The `Pileup`, or an error message.
    fn pileup(&mut self, chrom: &[u8], start: u32, variant: Variant) -> Result<Pileup<A, B, P>, String> {
        let case_pileup = try!(self.case_sample_mut().extract_observations(chrom, start, variant));
        let control_pileup = try!(self.control_sample_mut().extract_observations(chrom, start, variant));
        debug!("Obtained pileup (case: {} observations, control: {} observations).", case_pileup.len(), control_pileup.len());
        Ok(Pileup::new(
            case_pileup,
            control_pileup,
            variant
        ))
    }
}


/// Joint variant calling model, combining two latent variable models.
pub struct ContinuousVsDiscreteModel<P: priors::PairModel<ContinousAlleleFreqs, DiscreteAlleleFreqs>> {
    case_sample: Sample,
    control_sample: Sample,
    prior_model: P
}

impl<P: priors::PairModel<ContinousAlleleFreqs, DiscreteAlleleFreqs>> ContinuousVsDiscreteModel<P> {
    /// Create new `JointModel`.
    ///
    /// # Arguments
    ///
    /// * `case_sample` - case sample
    /// * `control_sample` - control sample
    pub fn new(case_sample: Sample, control_sample: Sample, prior_model: P) -> Self {
        ContinuousVsDiscreteModel {
            case_sample: case_sample,
            control_sample: control_sample,
            prior_model: prior_model
        }
    }
}


impl<P: Sync + priors::PairModel<ContinousAlleleFreqs, DiscreteAlleleFreqs>> JointModel<ContinousAlleleFreqs, DiscreteAlleleFreqs, P> for ContinuousVsDiscreteModel<P> {

    fn case_sample(&self) -> &Sample {
        &self.case_sample
    }

    fn control_sample(&self) -> &Sample {
        &self.control_sample
    }

    fn case_sample_mut(&mut self) -> &mut Sample {
        &mut self.case_sample
    }

    fn control_sample_mut(&mut self) -> &mut Sample {
        &mut self.control_sample
    }

    fn prior_model(&self) -> &P {
        &self.prior_model
    }

    fn joint_prob(
        &self,
        case_pileup: &[Observation],
        control_pileup: &[Observation],
        af_case: &ContinousAlleleFreqs,
        af_control: &DiscreteAlleleFreqs,
        variant: Variant
    ) -> LogProb {
        let case_likelihood = |af_case: AlleleFreq, af_control: AlleleFreq| {
            self.case_sample.likelihood_model().likelihood_pileup(case_pileup, *af_case, *af_control)
        };
        let control_likelihood = |af_control: AlleleFreq, _| {
            self.control_sample.likelihood_model().likelihood_pileup(control_pileup, *af_control, 0.0)
        };

        self.prior_model.joint_prob(af_case, af_control, &case_likelihood, &control_likelihood, variant)
    }
}


/// Pileup of observations associated with marginal probability.
pub struct Pileup<A, B, P> where
    A: AlleleFreqs,
    B: AlleleFreqs,
    P: priors::PairModel<A, B>
{
    case: Vec<Observation>,
    control: Vec<Observation>,
    // we use Cell for marginal prob to be able to mutate the field without having mutable access to the whole pileup
    marginal_prob: Cell<Option<LogProb>>,
    variant: Variant,
    a: PhantomData<A>,
    b: PhantomData<B>,
    p: PhantomData<P>
}


impl<A: AlleleFreqs, B: AlleleFreqs, P: priors::PairModel<A, B>> Pileup<A, B, P> {
    /// Create new pileup.
    fn new(case: Vec<Observation>, control: Vec<Observation>, variant: Variant) -> Self {
        Pileup {
            case: case,
            control: control,
            marginal_prob: Cell::new(None),
            variant: variant,
            a: PhantomData,
            b: PhantomData,
            p: PhantomData
        }
    }

    fn marginal_prob<M: JointModel<A, B, P>>(&self, model: &M) -> LogProb {
        if self.marginal_prob.get().is_none() {
            debug!("Calculating marginal probability.");

            self.marginal_prob.set(Some(model.marginal_prob(&self.case, &self.control, self.variant)));
        }

        self.marginal_prob.get().unwrap()
    }

    /// Calculate posterior probability of given allele frequencies.
    pub fn posterior_prob<M: JointModel<A, B, P>>(&self, model: &M, af_case: &A, af_control: &B) -> LogProb {
        debug!("Calculating posterior probability. for case={:?} and control={:?}", af_case, af_control);
        let p = model.joint_prob(&self.case, &self.control, af_case, af_control, self.variant);
        let marginal = self.marginal_prob(model);
        let prob = p - marginal;
        prob
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use likelihood::LatentVariableModel;
    use Sample;
    use InsertSize;
    use model::sample::Observation;
    use rust_htslib::bam;
    use bio::stats::LogProb;
    #[cfg(feature="flame_it")]
    use std::fs::File;
    #[cfg(feature="flame_it")]
    use flame;

    #[test]
    fn test_joint_prob() {
        let variant = Variant::Deletion(3);
        let insert_size = InsertSize{ mean: 250.0, sd: 50.0 };
        let prior_model = priors::TumorNormalModel::new(2, 30.0, 1.0, 1.0, 3e9 as u64, 1.0, 0.001);
        let case_sample = Sample::new(
            bam::IndexedReader::new(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            insert_size,
            LatentVariableModel::new(1.0)
        );
        let control_sample = Sample::new(
            bam::IndexedReader::new(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            true,
            insert_size,
            LatentVariableModel::new(1.0)
        );

        let model = ContinuousVsDiscreteModel::new(
            case_sample,
            control_sample,
            prior_model
        );

        let mut observations = Vec::new();
        for _ in 0..5 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_one(),
                prob_ref: LogProb::ln_zero(),
                prob_mismapped: LogProb::ln_one()
            });
        }

        let tumor_all = AlleleFreq(0.0)..AlleleFreq(1.0);
        let tumor_alt = AlleleFreq(0.001)..AlleleFreq(1.0);
        let tumor_ref = AlleleFreq(0.0)..AlleleFreq(0.001);
        let normal_alt = vec![AlleleFreq(0.5), AlleleFreq(1.0)];
        let normal_ref = vec![AlleleFreq(0.0)];

        // scenario 1: same pileup -> germline call
        let marginal_prob = model.marginal_prob(&observations, &observations, variant);
        let pileup = Pileup::new(observations.clone(), observations.clone(), variant);
        pileup.marginal_prob.set(Some(marginal_prob));
        // germline
        let p_germline = pileup.posterior_prob(&model, &tumor_all, &normal_alt);
        let p_somatic = pileup.posterior_prob(&model, &tumor_alt, &normal_ref);
        assert_relative_eq!(p_germline.exp(), 1.0);
        // somatic
        assert_relative_eq!(p_somatic.exp(), 0.0);
        assert!(p_germline.ln_add_exp(p_somatic).exp() <= 1.0);

        // scenario 2: empty control pileup -> somatic call
        let marginal_prob = model.marginal_prob(&observations, &[], variant);
        let pileup = Pileup::new(observations.clone(), vec![], variant);
        pileup.marginal_prob.set(Some(marginal_prob));
        let p_germline = pileup.posterior_prob(&model, &tumor_all, &normal_alt);
        let p_somatic = pileup.posterior_prob(&model, &tumor_alt, &normal_ref);
        let p_absent = pileup.posterior_prob(&model, &tumor_ref, &normal_ref);
        // without knowledge about germline, the germline prior dominates
        assert!(p_germline > p_somatic);
        assert!(p_germline > p_absent);
        // germline < somatic
        assert!(p_germline.ln_add_exp(p_somatic).exp() <= 1.0);

        // scenario 3: subclonal variant
        for _ in 0..50 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_zero(),
                prob_ref: LogProb::ln_one(),
                prob_mismapped: LogProb::ln_one()
            });
        }
        let marginal_prob = model.marginal_prob(&observations, &[], variant);
        let pileup = Pileup::new(observations.clone(), vec![], variant);
        pileup.marginal_prob.set(Some(marginal_prob));
        let p_germline = pileup.posterior_prob(&model, &tumor_all, &normal_alt);
        let p_somatic = pileup.posterior_prob(&model, &tumor_alt, &normal_ref);
        // somatic
        assert_relative_eq!(p_somatic.exp(), 0.9985, epsilon=0.01);
        assert!(p_germline.ln_add_exp(p_somatic).exp() <= 1.0);

        // scenario 4: absent variant
        observations.clear();
        for _ in 0..10 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_zero(),
                prob_ref: LogProb::ln_one(),
                prob_mismapped: LogProb::ln_one()
            });
        }
        println!("Absent call.");
        let marginal_prob = model.marginal_prob(&observations, &observations, variant);
        let pileup = Pileup::new(observations.clone(), observations.clone(), variant);
        pileup.marginal_prob.set(Some(marginal_prob));
        let p_absent = pileup.posterior_prob(&model, &tumor_ref, &normal_ref);
        let p_somatic = pileup.posterior_prob(&model, &tumor_alt, &normal_ref);
        let p_germline = pileup.posterior_prob(&model, &tumor_all, &normal_alt);

        // germline
        assert_relative_eq!(p_germline.exp(), 0.0, epsilon=0.01);
        // somatic
        assert_relative_eq!(p_somatic.exp(), 0.0, epsilon=0.01);
        assert_relative_eq!(p_absent.exp(), 1.0, epsilon=0.01);
    }

    #[cfg_attr(feature="flame_it", test)]
    #[cfg(feature="flame_it")]
    fn profile() {
        let variant = Variant::Deletion(3);
        let insert_size = InsertSize{ mean: 250.0, sd: 50.0 };
        let case_sample = Sample::new(
            bam::IndexedReader::new(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            insert_size,
            priors::TumorModel::new(2, 30.0, 1.0, 1.0, 3e9 as u64, 1.0, 0.001),
            LatentVariableModel::new(1.0)
        );
        let control_sample = Sample::new(
            bam::IndexedReader::new(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            true,
            insert_size,
            priors::InfiniteSitesNeutralVariationModel::new(2, 0.001),
            LatentVariableModel::new(1.0)
        );

        let model = ContinuousVsDiscreteModel::new(
            case_sample,
            control_sample
        );

        let mut observations = Vec::new();
        for _ in 0..1 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_one(),
                prob_ref: LogProb::ln_zero(),
                prob_mismapped: LogProb::ln_one()
            });
        }

        let tumor_all = 0.0..1.0;
        let tumor_alt = 0.001..1.0;
        let normal_alt = vec![0.5, 1.0];
        let normal_ref = vec![0.0];

        let marginal_prob = model.marginal_prob(&observations, &observations, variant);
    }
}

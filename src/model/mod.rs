use std::marker::PhantomData;

pub mod likelihood;
pub mod priors;
pub mod sample;

use bio::stats::LogProb;

use model::sample::{Sample, Observation};
use model::priors::{AlleleFreq, ContinousAlleleFreq, DiscreteAlleleFreq};


#[derive(Copy, Clone)]
pub enum Variant {
    Deletion(u32),
    Insertion(u32)
}


pub trait JointModel<A: AlleleFreq, B: AlleleFreq, P: priors::Model<A>, Q: priors::Model<B>> {
    /// Return case sample.
    fn case_sample(&self) -> &Sample<A, P>;

    /// Return control sample.
    fn control_sample(&self) -> &Sample<B, Q>;

    /// Return case sample.
    fn case_sample_mut(&mut self) -> &mut Sample<A, P>;

    /// Return control sample.
    fn control_sample_mut(&mut self) -> &mut Sample<B, Q>;

    /// Calculate joint probability of given pileups and allele frequencies.
    fn joint_prob(&self, case_pileup: &[Observation], control_pileup: &[Observation], af_case: &A, af_control: &B, variant: Variant) -> LogProb;

    /// Calculate marginal probability of given pileups.
    fn marginal_prob(&self, case_pileup: &[Observation], control_pileup: &[Observation], variant: Variant) -> LogProb {
        let af_case = self.case_sample().prior_model().allele_freqs();
        let af_control = self.control_sample().prior_model().allele_freqs();
        self.joint_prob(case_pileup, control_pileup, af_case, af_control, variant)
    }

    /// Calculate pileup and marginal probability for given variant.
    ///
    /// # Arguments
    ///
    /// * `chrom` - the chromosome of the variant
    /// * `start` - the starting position of the variant
    /// * `length` - the length of the variant
    /// * `is_del` - whether the variant is a deletion or an insertion
    ///
    /// # Returns
    /// The `Pileup`, or an error message.
    fn pileup(&mut self, chrom: &[u8], start: u32, variant: Variant) -> Result<Pileup<A, B, P, Q, Self>, String> where
        Self: Sized
    {
        let case_pileup = try!(self.case_sample_mut().extract_observations(chrom, start, variant));
        let control_pileup = try!(self.control_sample_mut().extract_observations(chrom, start, variant));
        Ok(Pileup::new(
            self,
            case_pileup,
            control_pileup,
            variant
        ))
    }
}


/// Joint variant calling model, combining two latent variable models.
pub struct ContinuousVsDiscreteModel<P: priors::Model<ContinousAlleleFreq>, Q: priors::Model<DiscreteAlleleFreq>> {
    case_sample: Sample<ContinousAlleleFreq, P>,
    control_sample: Sample<DiscreteAlleleFreq, Q>
}

impl<P: priors::Model<ContinousAlleleFreq>, Q: priors::Model<DiscreteAlleleFreq>> ContinuousVsDiscreteModel<P, Q> {
    /// Create new `JointModel`.
    ///
    /// # Arguments
    ///
    /// * `case_sample` - case sample
    /// * `control_sample` - control sample
    pub fn new(case_sample: Sample<ContinousAlleleFreq, P>, control_sample: Sample<DiscreteAlleleFreq, Q>) -> Self {
        ContinuousVsDiscreteModel {
            case_sample: case_sample,
            control_sample: control_sample
        }
    }
}


impl<P: priors::Model<ContinousAlleleFreq>, Q: priors::Model<DiscreteAlleleFreq>> JointModel<ContinousAlleleFreq, DiscreteAlleleFreq, P, Q> for ContinuousVsDiscreteModel<P, Q> {

    fn case_sample(&self) -> &Sample<ContinousAlleleFreq, P> {
        &self.case_sample
    }

    fn control_sample(&self) -> &Sample<DiscreteAlleleFreq, Q> {
        &self.control_sample
    }

    fn case_sample_mut(&mut self) -> &mut Sample<ContinousAlleleFreq, P> {
        &mut self.case_sample
    }

    fn control_sample_mut(&mut self) -> &mut Sample<DiscreteAlleleFreq, Q> {
        &mut self.control_sample
    }

    fn joint_prob(
        &self,
        case_pileup: &[Observation],
        control_pileup: &[Observation],
        af_case: &priors::ContinousAlleleFreq,
        af_control: &priors::DiscreteAlleleFreq,
        variant: Variant
    ) -> LogProb {
        let control_event_prob = |af_control| {
            let case_likelihood = |af| self.case_sample.likelihood_model().likelihood_pileup(case_pileup, af, af_control);
            let p_case = self.case_sample.prior_model().joint_prob(af_case, &case_likelihood, variant);
            let l_control = self.control_sample.likelihood_model().likelihood_pileup(control_pileup, af_control, 0.0);

            p_case + l_control
        };

        let prob = self.control_sample.prior_model().joint_prob(af_control, &control_event_prob, variant);

        prob
    }
}


/// Pileup of observations associated with marginal probability.
pub struct Pileup<'a, A, B, P, Q, M> where
    A: AlleleFreq,
    B: AlleleFreq,
    P: priors::Model<A>,
    Q: priors::Model<B>,
    M: 'a + JointModel<A, B, P, Q>
{
    model: &'a M,
    case: Vec<Observation>,
    control: Vec<Observation>,
    marginal_prob: Option<LogProb>,
    variant: Variant,
    a: PhantomData<A>,
    b: PhantomData<B>,
    p: PhantomData<P>,
    q: PhantomData<Q>
}


impl<'a, A: AlleleFreq, B: AlleleFreq, P: priors::Model<A>, Q: priors::Model<B>, M: JointModel<A, B, P, Q>> Pileup<'a, A, B, P, Q, M> {
    /// Create new pileup.
    fn new(model: &'a M, case: Vec<Observation>, control: Vec<Observation>, variant: Variant) -> Self {
        Pileup {
            model: model,
            case: case,
            control: control,
            marginal_prob: None,
            variant: variant,
            a: PhantomData,
            b: PhantomData,
            p: PhantomData,
            q: PhantomData
        }
    }

    fn marginal_prob(&mut self) -> LogProb {
        if self.marginal_prob.is_none() {
            debug!("Calculating marginal probability.");
            self.marginal_prob = Some(self.model.marginal_prob(&self.case, &self.control, self.variant));
        }

        self.marginal_prob.unwrap()
    }

    /// Calculate posterior probability of given allele frequencies.
    pub fn posterior_prob(&mut self, af_case: &A, af_control: &B) -> LogProb {
        debug!("Calculating posterior probability");
        let prob = self.model.joint_prob(&self.case, &self.control, af_case, af_control, self.variant) - self.marginal_prob();

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

    #[test]
    fn test_joint_prob() {
        let variant = Variant::Deletion(3);
        let insert_size = InsertSize{ mean: 250.0, sd: 50.0 };
        let case_sample = Sample::new(
            bam::IndexedReader::new(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            insert_size,
            priors::TumorModel::new(2, 30.0, 1.0, 1.0, 3e9 as u64, 1.0, 0.001),
            LatentVariableModel::new(1.0)
        );
        let control_sample = Sample::new(
            bam::IndexedReader::new(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            insert_size,
            priors::InfiniteSitesNeutralVariationModel::new(2, 0.001),
            LatentVariableModel::new(1.0)
        );

        let model = ContinuousVsDiscreteModel::new(
            case_sample,
            control_sample
        );

        let mut observations = Vec::new();
        for _ in 0..5 {
            observations.push(Observation{
                prob_mapping: 1.0f64.ln(),
                prob_alt: 1.0f64.ln(),
                prob_ref: 0.0f64.ln(),
                prob_mismapped: 1.0f64.ln()
            });
        }

        let tumor_all = 0.0..1.0;
        let tumor_alt = 0.001..1.0;
        let normal_alt = vec![0.5, 1.0];
        let normal_ref = vec![0.0];

        // scenario 1: same pileup -> germline call
        let marginal_prob = model.marginal_prob(&observations, &observations, variant);
        let pileup = Pileup::new(&model, observations.clone(), observations.clone(), marginal_prob, variant);
        // germline
        assert_relative_eq!(pileup.posterior_prob(&tumor_all, &normal_alt).exp(), 1.0);
        // somatic
        assert_relative_eq!(pileup.posterior_prob(&tumor_alt, &normal_ref).exp(), 0.0);

        // scenario 2: empty control pileup -> somatic call
        let marginal_prob = model.marginal_prob(&observations, &[], variant);
        let pileup = Pileup::new(&model, observations.clone(), vec![], marginal_prob, variant);
        // somatic close to prior for ref in ńormal
        assert_relative_eq!(pileup.posterior_prob(&tumor_alt, &normal_ref).exp(), 0.9985, epsilon=0.01);
        // germline < somatic
        assert!(pileup.posterior_prob(&tumor_all, &normal_alt).exp() < pileup.posterior_prob(&tumor_alt, &normal_ref).exp());

        // scenario 3: subclonal variant
        for _ in 0..50 {
            observations.push(Observation{
                prob_mapping: 1.0f64.ln(),
                prob_alt: 0.0f64.ln(),
                prob_ref: 1.0f64.ln(),
                prob_mismapped: 1.0f64.ln()
            });
        }
        let marginal_prob = model.marginal_prob(&observations, &[], variant);
        let pileup = Pileup::new(&model, observations.clone(), vec![], marginal_prob, variant);
        // somatic
        assert_relative_eq!(pileup.posterior_prob(&tumor_alt, &normal_ref).exp(), 0.9985, epsilon=0.01);

        // scenario 4: absent variant
        observations.clear();
        for _ in 0..10 {
            observations.push(Observation{
                prob_mapping: 1.0f64.ln(),
                prob_alt: 0.0f64.ln(),
                prob_ref: 1.0f64.ln(),
                prob_mismapped: 1.0f64.ln()
            });
        }
        let marginal_prob = model.marginal_prob(&observations, &observations, variant);
        let pileup = Pileup::new(&model, observations.clone(), observations.clone(), marginal_prob, variant);
        // germline
        assert_relative_eq!(pileup.posterior_prob(&tumor_all, &normal_alt).exp(), 0.0, epsilon=0.01);
        // somatic
        assert_relative_eq!(pileup.posterior_prob(&tumor_alt, &normal_ref).exp(), 0.0, epsilon=0.01);
    }
}

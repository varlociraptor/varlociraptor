use std::ops::Range;

pub mod likelihood;
pub mod priors;
pub mod sample;

use bio::stats::{LogProb, logprobs};

use model::sample::{Sample, Observation};


#[derive(Copy, Clone)]
pub enum Variant {
    Deletion(u32),
    Insertion(u32)
}


/// Joint variant calling model, combining two latent variable models.
pub struct JointModel<P: priors::Model, Q: priors::Model> {
    case_sample: Sample<P>,
    control_sample: Sample<Q>
}


impl<P: priors::ContinuousModel, Q: priors::DiscreteModel> JointModel<P, Q> {

    /// Create new `JointModel`.
    ///
    /// # Arguments
    ///
    /// * `case_sample` - case sample
    /// * `control_sample` - control sample
    /// * `grid_points` - number of grid points to use for trapezoidal integration (e.g. 200)
    pub fn new(case_sample: Sample<P>, control_sample: Sample<Q>) -> Self {
        JointModel {
            case_sample: case_sample,
            control_sample: control_sample
        }
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
    pub fn pileup(&mut self, chrom: &[u8], start: u32, variant: Variant) -> Result<Pileup<P, Q>, String> {
        let case_pileup = try!(self.case_sample.extract_observations(chrom, start, variant));
        let control_pileup = try!(self.control_sample.extract_observations(chrom, start, variant));
        let marginal_prob = self.marginal_prob(&case_pileup, &control_pileup);
        Ok(Pileup::new(
            self,
            case_pileup,
            control_pileup,
            marginal_prob
        ))
    }

    fn joint_prob(&self, case_pileup: &[Observation], control_pileup: &[Observation], af_case: &Range<f64>, af_control: f64) -> LogProb {
        let case_likelihood = |af| self.case_sample.likelihood_model().likelihood_pileup(case_pileup, af, af_control);

        let prob = self.control_sample.prior_model().prior_prob(af_control) +
                   self.control_sample.likelihood_model().likelihood_pileup(control_pileup, af_control, 0.0) +
                   self.case_sample.prior_model().integrate(af_case, &case_likelihood);
        prob
    }

    fn marginal_prob(&self, case_pileup: &[Observation], control_pileup: &[Observation]) -> LogProb {
        let mut summands = Vec::with_capacity(3);
        for &af_control in [0.0, 0.5, 1.0].iter() {
            summands.push(self.joint_prob(case_pileup, control_pileup, &(0.0..1.0), af_control));
        }

        logprobs::sum(&summands)
    }
}


/// Pileup of observations associated with marginal probability.
pub struct Pileup<'a, P: 'a + priors::Model, Q: 'a + priors::Model> {
    model: &'a JointModel<P, Q>,
    case: Vec<Observation>,
    control: Vec<Observation>,
    marginal_prob: LogProb
}


impl<'a, P: priors::Model, Q: priors::Model> Pileup<'a, P, Q> {
    /// Create new pileup.
    fn new(model: &'a JointModel<P, Q>, case: Vec<Observation>, control: Vec<Observation>, marginal_prob: LogProb) -> Self {
        Pileup {
            model: model,
            case: case,
            control: control,
            marginal_prob: marginal_prob
        }
    }
}

impl<'a, P: priors::ContinuousModel, Q: priors::DiscreteModel> Pileup<'a, P, Q> {
    /// Calculate posterior probability of given allele frequencies.
    pub fn posterior_prob(&self, af_case: &Range<f64>, af_control: &[f64]) -> LogProb {
        let mut summands = Vec::with_capacity(af_control.len());
        for &af_control in af_control.iter() {
            summands.push(self.model.joint_prob(&self.case, &self.control, af_case, af_control));
        }
        logprobs::sum(&summands) - self.marginal_prob
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use LatentVariableModel;
    use Sample;
    use InsertSize;
    use model::sample::Observation;
    use rust_htslib::bam;

    #[test]
    fn test_joint_prob() {
        let insert_size = InsertSize{ mean: 250.0, sd: 50.0 };
        let case_sample = Sample::new(
            bam::IndexedReader::new(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            insert_size,
            priors::TumorModel::new(2, 30.0, 3e9 as u64, 1.0, 0.001)
        );
        let control_sample = Sample::new(
            bam::IndexedReader::new(&"tests/test.bam").expect("Error reading BAM."),
            5000,
            insert_size,
            priors::InfiniteSitesNeutralVariationModel::new(2, 0.001)
        );

        let model = JointModel::new(
            LatentVariableModel::new(1.0),
            LatentVariableModel::new(1.0),
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
        let normal_alt = [0.5, 1.0];
        let normal_ref = [0.0];

        // scenario 1: same pileup -> germline call
        let marginal_prob = model.marginal_prob(&observations, &observations);
        let pileup = Pileup::new(&model, observations.clone(), observations.clone(), marginal_prob);
        // germline
        assert_relative_eq!(pileup.posterior_prob(&tumor_all, &normal_alt).exp(), 1.0);
        // somatic
        assert_relative_eq!(pileup.posterior_prob(&tumor_alt, &normal_ref).exp(), 0.0);

        // scenario 2: empty control pileup -> somatic call
        let marginal_prob = model.marginal_prob(&observations, &[]);
        let pileup = Pileup::new(&model, observations.clone(), vec![], marginal_prob);
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
        let marginal_prob = model.marginal_prob(&observations, &[]);
        let pileup = Pileup::new(&model, observations.clone(), vec![], marginal_prob);
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
        let marginal_prob = model.marginal_prob(&observations, &observations);
        let pileup = Pileup::new(&model, observations.clone(), observations.clone(), marginal_prob);
        // germline
        assert_relative_eq!(pileup.posterior_prob(&tumor_all, &normal_alt).exp(), 0.0, epsilon=0.01);
        // somatic
        assert_relative_eq!(pileup.posterior_prob(&tumor_alt, &normal_ref).exp(), 0.0, epsilon=0.01);
    }
}

use std::ops::Range;

pub mod likelihood;
pub mod priors;
pub mod sample;

use bio::stats::{LogProb, logprobs};

use model::likelihood::LatentVariableModel;
use model::sample::{Sample, Observation};


#[derive(Copy, Clone)]
pub enum Variant {
    Deletion(u32),
    Insertion(u32)
}


/// Joint variant calling model, combining two latent variable models.
pub struct JointModel<P: priors::Model, Q: priors::Model> {
    case_model: LatentVariableModel,
    control_model: LatentVariableModel,
    case_sample: Sample<P>,
    control_sample: Sample<Q>,
    grid_points: usize
}


impl<P: priors::ContinuousModel, Q: priors::DiscreteModel> JointModel<P, Q> {

    /// Create new `JointModel`.
    ///
    /// # Arguments
    ///
    /// * `case_model` - model for the case sample
    /// * `control_model` - model for the control sample
    /// * `case_sample` - case sample
    /// * `control_sample` - control sample
    /// * `grid_points` - number of grid points to use for trapezoidal integration (e.g. 200)
    pub fn new(case_model: LatentVariableModel, control_model: LatentVariableModel, case_sample: Sample<P>, control_sample: Sample<Q>, grid_points: usize) -> Self {
        JointModel {
            case_model: case_model,
            control_model: control_model,
            case_sample: case_sample,
            control_sample: control_sample,
            grid_points: grid_points
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
        let case_density = |af| {
            self.case_model.likelihood_pileup(case_pileup, af, af_control) +
            self.case_sample.prior_prob(af)
        };

        let prob = self.control_sample.prior_prob(af_control) +
                   self.control_model.likelihood_pileup(control_pileup, af_control, 0.0) +
                   logprobs::integrate(case_density, af_case.start, af_case.end, self.grid_points);
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

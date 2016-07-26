use std::ops::Range;

pub mod likelihood;
pub mod priors;
pub mod sample;

use bio::stats::{LogProb, Prob, logprobs};
use rgsl;
use rgsl::integration;

use model::likelihood::LatentVariableModel;
use model::sample::{Sample, Observation};


/// Joint variant calling model, combining two latent variable models.
pub struct JointModel<P: priors::Model, Q: priors::Model> {
    case_model: LatentVariableModel,
    control_model: LatentVariableModel,
    case_sample: Sample<P>,
    control_sample: Sample<Q>
}


impl<P: priors::ContinuousModel, Q: priors::DiscreteModel> JointModel<P, Q> {

    /// Create new `JointModel`.
    pub fn new(case_model: LatentVariableModel, control_model: LatentVariableModel, case_sample: Sample<P>, control_sample: Sample<Q>) -> Self {
        JointModel {
            case_model: case_model,
            control_model: control_model,
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
    pub fn pileup(&mut self, chrom: &[u8], start: u32, length: u32, is_del: bool) -> Result<Pileup<P, Q>, String> {
        let case_pileup = try!(self.case_sample.extract_observations(chrom, start, length, is_del));
        let control_pileup = try!(self.control_sample.extract_observations(chrom, start, length, is_del));
        let marginal_prob = try!(self.marginal_prob(&case_pileup, &control_pileup));
        Ok(Pileup::new(
            self,
            case_pileup,
            control_pileup,
            marginal_prob
        ))
    }

    fn prior_prob(&self, af_case: f64, af_control: f64) -> LogProb {
        self.case_sample.prior_prob(af_case) + self.control_sample.prior_prob(af_control)
    }

    fn joint_prob(&self, case_pileup: &[Observation], control_pileup: &[Observation], af_case: &Range<f64>, af_control: f64) -> Result<LogProb, String> {
        fn f<P: priors::ContinuousModel, Q: priors::DiscreteModel>(af_case: f64, params: &mut (&JointModel<P, Q>, &[Observation], f64, LogProb)) -> Prob {
            let (_self, case_pileup, af_control, lh_control) = *params;
            (_self.case_model.likelihood_pileup(case_pileup, af_case, af_control) +
            lh_control + _self.prior_prob(af_case, af_control)).exp()
        }

        let lh_control = self.control_model.likelihood_pileup(control_pileup, af_control, 0.0);
        let mut prob = 0.0f64;
        let mut abs_err = 0.0;
        let mut n_eval = 0;
        if let rgsl::Value::Success = integration::qng(f, &mut (&self, case_pileup, af_control, lh_control), af_case.start, af_case.end, 1e-8, 1e-8, &mut prob, &mut abs_err, &mut n_eval) {
            Ok(prob)
        } else {
            return Err("Error calculating integral.".to_owned())
        }
    }

    fn marginal_prob(&self, case_pileup: &[Observation], control_pileup: &[Observation]) -> Result<LogProb, String> {
        let mut summands = Vec::with_capacity(3);
        for &af_control in [0.0, 0.5, 1.0].iter() {
            summands.push(try!(self.joint_prob(case_pileup, control_pileup, &(0.0..1.0), af_control)).ln());
        }

        Ok(logprobs::sum(&summands))
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
    pub fn posterior_prob(&self, af_case: &Range<f64>, af_control: &[f64]) -> Result<LogProb, String> {
        let mut summands = Vec::with_capacity(af_control.len());
        for &af_control in af_control.iter() {
            summands.push(try!(self.model.joint_prob(&self.case, &self.control, af_case, af_control)));
        }
        Ok(logprobs::sum(&summands) - self.marginal_prob)
    }
}

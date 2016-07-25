use std::ops::Range;

pub mod likelihood;
pub mod observations;
pub mod priors;

use bio::stats::{LogProb, Prob, logprobs};
use rgsl;
use rgsl::integration;

use model::likelihood::LatentVariableModel;
use model::observations::Observation;


#[derive(Copy, Clone, Debug)]
pub struct InsertSize {
    mean: f64,
    sd: f64
}


pub struct JointModel<P: priors::Model, Q: priors::Model> {
    case_model: LatentVariableModel,
    control_model: LatentVariableModel,
    case_prior_model: P,
    control_prior_model: Q
}


impl<P: priors::Model, Q: priors::Model> JointModel<P, Q> {

    pub fn new(case_model: LatentVariableModel, control_model: LatentVariableModel, case_prior_model: P, control_prior_model: Q) -> Self {
        JointModel {
            case_model: case_model,
            control_model: control_model,
            case_prior_model: case_prior_model,
            control_prior_model: control_prior_model
        }
    }

    fn prior_prob(&self, af_case: f64, af_control: f64) -> LogProb {
        self.case_prior_model.prior_prob(af_case) + self.control_prior_model.prior_prob(af_control)
    }

    fn joint_prob(&self, case_pileup: &[Observation], control_pileup: &[Observation], af_case: &Range<f64>, af_control: f64) -> Result<LogProb, String> {
        fn f<P: priors::Model, Q: priors::Model>(af_case: f64, params: &mut (&JointModel<P, Q>, &[Observation], f64, LogProb)) -> Prob {
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

    pub fn marginal_prob(&self, case_pileup: &[Observation], control_pileup: &[Observation]) -> Result<LogProb, String> {
        let mut summands = Vec::with_capacity(3);
        for &af_control in [0.0, 0.5, 1.0].iter() {
            summands.push(try!(self.joint_prob(case_pileup, control_pileup, &(0.0..1.0), af_control)).ln());
        }

        Ok(logprobs::sum(&summands))
    }

    pub fn posterior_prob(&self, case_pileup: &[Observation], control_pileup: &[Observation], af_case: &Range<f64>, af_control: &[f64], marginal_prob: LogProb) -> Result<LogProb, String> {
        let mut summands = Vec::with_capacity(af_control.len());
        for &af_control in af_control.iter() {
            summands.push(try!(self.joint_prob(case_pileup, control_pileup, af_case, af_control)));
        }
        Ok(logprobs::sum(&summands) - marginal_prob)
    }
}

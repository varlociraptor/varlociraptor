use std::ops::Range;

pub mod likelihood;
pub mod observations;
pub mod priors;
pub mod sample;

use bio::stats::{LogProb, Prob, logprobs};
use rgsl;
use rgsl::integration;

use model::likelihood::LatentVariableModel;
use model::sample::{Sample, Observation};


#[derive(Copy, Clone, Debug)]
pub struct InsertSize {
    mean: f64,
    sd: f64
}


pub struct JointModel<A: Sample, B: Sample> {
    case_model: LatentVariableModel,
    control_model: LatentVariableModel,
    case_sample: A,
    control_sample: B
}


impl<A: Sample, B: Sample> JointModel<A, B> {

    pub fn new(case_model: LatentVariableModel, control_model: LatentVariableModel, case_sample: A, control_sample: B) -> Self {
        JointModel {
            case_model: case_model,
            control_model: control_model,
            case_sample: case_sample,
            control_sample: control_sample
        }
    }

    pub fn pileup(&mut self, chrom: &[u8], start: u32, length: u32, is_del: bool) -> Result<Pileup<A, B>, String> {
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
        fn f<A: Sample, B: Sample>(af_case: f64, params: &mut (&JointModel<A, B>, &[Observation], f64, LogProb)) -> Prob {
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


pub struct Pileup<'a, A: 'a + Sample, B: 'a + Sample> {
    model: &'a JointModel<A, B>,
    case: Vec<Observation>,
    control: Vec<Observation>,
    marginal_prob: LogProb
}


impl<'a, A: Sample, B: Sample> Pileup<'a, A, B> {
    pub fn new(model: &'a JointModel<A, B>, case: Vec<Observation>, control: Vec<Observation>, marginal_prob: LogProb) -> Self {
        Pileup {
            model: model,
            case: case,
            control: control,
            marginal_prob: marginal_prob
        }
    }

    pub fn posterior_prob(&self, af_case: &Range<f64>, af_control: &[f64]) -> Result<LogProb, String> {
        let mut summands = Vec::with_capacity(af_control.len());
        for &af_control in af_control.iter() {
            summands.push(try!(self.model.joint_prob(&self.case, &self.control, af_case, af_control)));
        }
        Ok(logprobs::sum(&summands) - self.marginal_prob)
    }
}

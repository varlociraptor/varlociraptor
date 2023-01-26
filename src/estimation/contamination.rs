use std::{
    collections::BTreeMap,
    convert::TryFrom,
    fs::File,
    io::{stdout, Write},
    path::PathBuf,
};

use anyhow::Result;
use bio::stats::{bayesian, LogProb, Prob};
use csv::WriterBuilder;
use itertools::Itertools;
use itertools_num::linspace;
use rgsl::randist::binomial::binomial_pdf;
use serde_json;
use serde_json::json;

use crate::{
    calling::variants::{
        calling::{call_generic, CallProcessor, Caller, CandidateFilter, WorkItem},
        Call,
    },
    grammar,
    utils::aux_info::AuxInfoCollector,
    utils::PathMap,
    variants::model::AlleleFreq,
};

#[derive(Hash, Eq, PartialEq, Clone, Debug)]
struct Event {
    contamination: AlleleFreq,
    expected_max_somatic_vaf: AlleleFreq,
}

#[derive(Clone, Debug)]
struct VariantObservation {
    prob_denovo: LogProb,
    vaf_dist: BTreeMap<AlleleFreq, LogProb>,
    max_posterior_vaf: AlleleFreq,
}

impl VariantObservation {
    fn new(call: &Call, sample_names: &grammar::SampleInfo<String>) -> Option<Self> {
        let sample_info = call.variant().as_ref().unwrap().sample_info()
            [sample_names.iter().position(|s| *s == "sample").unwrap()]
        .as_ref()
        .unwrap();
        let prob_denovo = *call
            .variant()
            .as_ref()
            .unwrap()
            .event_probs()
            .as_ref()
            .unwrap()
            .get("denovo")
            .unwrap();

        if sample_info.vaf_dist().is_none() || prob_denovo.exp() < 0.95 {
            // no denovo variant, skip
            return None;
        }

        let mut vaf_dist = sample_info
            .vaf_dist()
            .as_ref()
            .unwrap()
            .iter()
            .map(|(vaf, density)| (*vaf, *density))
            .collect();

        Some(VariantObservation {
            vaf_dist,
            prob_denovo,
            max_posterior_vaf: sample_info.allelefreq_estimate(),
        })
    }

    fn pdf(&self, vaf: AlleleFreq) -> LogProb {
        let supremum = self.vaf_dist.range(vaf..).next();
        if let Some((supremum, density)) = supremum {
            if *supremum == vaf {
                // case 1: exact vaf found
                return *density;
            }
        }

        let infimum = self.vaf_dist.range(..vaf).last();

        match (infimum, supremum) {
            (Some((infimum, inf_density)), Some((supremum, sup_density))) => {
                // case 2: interpolate
                inf_density.ln_add_exp(LogProb(
                    ((sup_density.exp() - inf_density.exp()) / (**supremum - **infimum)).ln()
                        + (vaf - *infimum).ln(),
                ))
            }
            (Some((infimum, count)), None) => {
                // case 3: right of highest value, return zero
                LogProb::ln_zero()
            }
            (None, Some(_)) => {
                // case 4: left of smallest value, return zero
                LogProb::ln_zero()
            }
            (None, None) => {
                // case 5: empty dist
                LogProb::ln_zero()
            }
        }
    }
}

#[derive(new, Clone, Debug)]
struct Prior {
    prior_estimate: Option<PriorEstimate>,
}

impl Prior {
    fn as_json(&self) -> serde_json::Value {
        serde_json::Value::Array(
            linspace(0.0, 1.0, 101)
                .map(|contamination| {
                    json!({
                        "purity": 1.0 - contamination,
                        "density": self.prob(AlleleFreq(contamination)).exp(),
                        "category": "prior",
                    })
                })
                .collect_vec(),
        )
    }

    fn prob(&self, contamination: AlleleFreq) -> LogProb {
        if let Some(ref prior_estimate) = self.prior_estimate {
            let n = prior_estimate.n_observed_cells;
            let k = (prior_estimate.contamination * n as f64).round() as u32;
            let p = *contamination;
            LogProb::from(Prob(binomial_pdf(k, p, n)))
        } else {
            LogProb::ln_one()
        }
    }
}

impl bio::stats::bayesian::model::Prior for Prior {
    type Event = Event;

    fn compute(&self, event: &Self::Event) -> LogProb {
        self.prob(event.contamination)
    }
}

#[derive(new)]
struct Likelihood {
    vaf_dist: VAFDist,
}

impl bio::stats::bayesian::model::Likelihood for Likelihood {
    type Event = Event;
    type Data = Vec<VariantObservation>;

    fn compute(&self, event: &Self::Event, data: &Self::Data, payload: &mut ()) -> LogProb {
        let purity = AlleleFreq(1.0) - event.contamination;
        data.iter()
            .map(|obs| {
                if *purity == 0.0 {
                    // there cannot be any denovo somatic mutation
                    return obs.prob_denovo.ln_one_minus_exp();
                }

                let expected_vaf = self.vaf_dist.get_expected_vaf(
                    purity,
                    obs.max_posterior_vaf,
                    event.expected_max_somatic_vaf,
                );
                obs.pdf(expected_vaf)
            })
            .sum()
    }
}

#[derive(new)]
struct Posterior;

impl bio::stats::bayesian::model::Posterior for Posterior {
    type Event = Event;
    type Data = Vec<VariantObservation>;
    type BaseEvent = Event;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        joint_prob(event, data)
    }
}

#[derive(new)]
struct Marginal;

impl bio::stats::bayesian::model::Marginal for Marginal {
    type Event = Event;
    type Data = Vec<VariantObservation>;
    type BaseEvent = Event;

    fn compute<F: FnMut(&Self::Event, &Self::Data) -> LogProb>(
        &self,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        LogProb::ln_sum_exp(
            &[0.25, 0.5, 0.75, 1.0]
                .iter()
                .map(|expected_max_somatic_vaf| {
                    let density = |_, contamination| {
                        let event = Event {
                            contamination: AlleleFreq(contamination),
                            expected_max_somatic_vaf: AlleleFreq(*expected_max_somatic_vaf),
                        };
                        joint_prob(&event, data)
                    };
                    LogProb::ln_simpsons_integrate_exp(density, 0.0, 1.0, 101)
                })
                .collect_vec(),
        )
    }
}

#[derive(Clone, Debug)]
struct VAFDist {
    histogram: BTreeMap<AlleleFreq, usize>,
    total: usize,
    max_vaf: AlleleFreq,
}

impl VAFDist {
    fn new(variant_observations: &[VariantObservation]) -> Self {
        let mut histogram = BTreeMap::default();
        let mut total = 0;
        let mut max_vaf = AlleleFreq(0.0);
        for obs in variant_observations {
            let bin = AlleleFreq((*obs.max_posterior_vaf * 100.0).floor() / 100.0);
            *histogram.entry(bin).or_insert(0) += 1;
            total += 1;
            if obs.max_posterior_vaf > max_vaf {
                max_vaf = obs.max_posterior_vaf;
            }
        }
        VAFDist {
            histogram,
            total,
            max_vaf,
        }
    }

    fn quantiles(&self) -> Vec<AlleleFreq> {
        vec![self.max_vaf]
    }

    fn get_expected_vaf(
        &self,
        purity: AlleleFreq,
        max_posterior_vaf: AlleleFreq,
        expected_max_somatic_vaf: AlleleFreq,
    ) -> AlleleFreq {
        let quantile = max_posterior_vaf / self.max_vaf;
        expected_max_somatic_vaf * purity * quantile
    }

    fn hist_as_json(&self) -> serde_json::Value {
        serde_json::Value::Array(
            self.histogram
                .iter()
                .map(|(vaf, count)| {
                    json!({
                        "vaf": *vaf,
                        "count": *count
                    })
                })
                .collect_vec(),
        )
    }

    fn quantiles_as_json(&self) -> serde_json::Value {
        serde_json::Value::Array(
            self.quantiles()
                .iter()
                .map(|vaf| json!({"vaf": *vaf}))
                .collect_vec(),
        )
    }
}

#[derive(Debug, Clone, Copy, new)]
pub struct PriorEstimate {
    contamination: AlleleFreq,
    n_observed_cells: u32,
}

#[derive(Debug, Clone)]
pub(crate) struct ContaminationEstimator {
    prior_estimate: Option<PriorEstimate>,
    variant_observations: Vec<VariantObservation>,
    output: Option<PathBuf>,
    output_plot: Option<PathBuf>,
}

impl ContaminationEstimator {
    pub(crate) fn new(
        output: Option<PathBuf>,
        output_plot: Option<PathBuf>,
        prior_estimate: Option<PriorEstimate>,
    ) -> Self {
        ContaminationEstimator {
            variant_observations: Vec::new(),
            output,
            output_plot,
            prior_estimate,
        }
    }

    fn calc_posterior<W: Write>(&self, mut writer: csv::Writer<W>) -> Result<()> {
        let vaf_dist = VAFDist::new(&self.variant_observations);

        let prior = Prior::new(self.prior_estimate);
        let likelihood = Likelihood::new(vaf_dist.clone());
        let model = bio::stats::bayesian::Model::new(likelihood, prior.clone(), Posterior::new());

        let model_instance =
            model.compute_from_marginal(&Marginal::new(), &self.variant_observations);

        if let Some(ref outpath) = self.output_plot {
            let mut spec = serde_json::from_str(include_str!(
                "../../templates/plots/contamination_estimation.json"
            ))?;
            let mut densities = prior.as_json();
            densities
                .as_array_mut()
                .unwrap()
                .extend(model_instance.event_posteriors().map(|(event, density)| {
                json!({
                    "purity": 1.0 - *event.contamination,
                    "density": *Prob::from(density),
                    "category": &format!("posterior, max VAF={}", *event.expected_max_somatic_vaf)
                })
            }));
            if let serde_json::Value::Object(ref mut spec) = spec {
                spec["datasets"]["empirical_vaf_dist"] = vaf_dist.hist_as_json();
                spec["datasets"]["vaf_dist_quantiles"] = vaf_dist.quantiles_as_json();
                spec["datasets"]["densities"] = densities;

                let mut outfile = File::create(outpath)?;
                serde_json::to_writer_pretty(outfile, &spec)?;
            } else {
                unreachable!();
            }
        }

        // write into table
        writer.write_record(&["maximum somatic VAF", "contamination", "posterior density"])?;
        for (event, density) in model_instance.event_posteriors() {
            writer.write_record(&[
                format!("{}", *event.expected_max_somatic_vaf),
                format!("{}", *event.contamination),
                format!("{}", *Prob::from(density)),
            ])?;
        }
        Ok(())
    }
}

impl CallProcessor for ContaminationEstimator {
    fn setup<Pr: bayesian::model::Prior, CF: CandidateFilter>(
        &mut self,
        caller: &Caller<Pr, Self, CF>,
    ) -> Result<Option<AuxInfoCollector>> {
        Ok(None)
    }

    fn process_call(
        &mut self,
        call: Call,
        sample_names: &grammar::SampleInfo<String>,
    ) -> Result<()> {
        let obs = VariantObservation::new(&call, sample_names);
        if let Some(obs) = obs {
            self.variant_observations.push(obs);
        }

        Ok(())
    }

    fn finalize(&mut self) -> Result<()> {
        if let Some(ref path) = self.output {
            self.calc_posterior(WriterBuilder::new().delimiter(b'\t').from_path(path)?)
        } else {
            self.calc_posterior(WriterBuilder::new().delimiter(b'\t').from_writer(stdout()))
        }
    }
}

#[derive(new)]
struct ContaminationCandidateFilter;

impl CandidateFilter for ContaminationCandidateFilter {
    fn filter(&self, work_item: &WorkItem, sample_names: &grammar::SampleInfo<String>) -> bool {
        // only keep variants where all reads of the contaminant are REF and that are SNVs
        let contaminant_pileup = &work_item.pileups().as_ref().unwrap()[sample_names
            .iter()
            .position(|s| *s == "contaminant")
            .unwrap()];
        let sample_pileup = &work_item.pileups().as_ref().unwrap()
            [sample_names.iter().position(|s| *s == "sample").unwrap()];
        work_item.snv().is_some()
            && contaminant_pileup.read_observations().len() >= 10
            && contaminant_pileup
                .read_observations()
                .iter()
                .all(|obs| obs.is_ref_support())
            && sample_pileup.read_observations().len() >= 10
            && sample_pileup
                .read_observations()
                .iter()
                .any(|obs| obs.is_strong_alt_support())
    }
}

pub(crate) fn estimate_contamination(
    sample: PathBuf,
    contaminant: PathBuf,
    output: Option<PathBuf>,
    output_plot: Option<PathBuf>,
    prior_estimate: Option<PriorEstimate>,
) -> Result<()> {
    let scenario = grammar::Scenario::try_from(
        r#"
    samples:
      sample:
        resolution: 0.01
        universe: "[0.0,1.0]"
      contaminant:
        resolution: 0.01
        universe: "[0.0,1.0]"
    events:
      denovo:  "sample:]0.0,1.0] & contaminant:0.0"
      other: "sample:[0.0,1.0] & contaminant:]0.0,1.0]"
    "#,
    )
    .unwrap();

    let mut observations = PathMap::default();
    observations.insert("sample".to_owned(), sample);
    observations.insert("contaminant".to_owned(), contaminant);

    call_generic(
        scenario,
        observations,
        false,
        false,
        false,
        false,
        false,
        false,
        None,
        false,
        ContaminationEstimator::new(output, output_plot, prior_estimate),
        ContaminationCandidateFilter::new(),
        Vec::new(),
    )
}

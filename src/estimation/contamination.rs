use std::{
    convert::TryFrom,
    io::{stdout, Write},
    path::PathBuf,
    collections::BTreeMap,
    fs::File,
};

use anyhow::Result;
use bio::stats::{bayesian, LogProb, Prob};
use csv::WriterBuilder;
use itertools::Itertools;
use itertools_num::linspace;
use ordered_float::NotNan;
use statrs::distribution::Empirical;
use serde_json::json;
use serde_json;

use crate::{
    calling::variants::{
        calling::{call_generic, CallProcessor, Caller, CandidateFilter, WorkItem},
        Call,
    },
    grammar,
    utils::PathMap,
    variants::model::AlleleFreq,
};

struct VariantObservation {
    prob_denovo: LogProb,
    vaf_dist: BTreeMap<AlleleFreq, LogProb>,
    max_posterior_vaf: AlleleFreq,
}

impl VariantObservation {
    fn new(call: &Call, sample_names: &grammar::SampleInfo<String>) -> Option<Self> {
        let sample_info = call.variant().as_ref().unwrap().sample_info()[sample_names.iter().position(|s| *s == "sample").unwrap()]
            .as_ref()
            .unwrap();
        let prob_denovo = *call.variant().as_ref().unwrap().event_probs().as_ref().unwrap().get("denovo").unwrap();

        if sample_info.vaf_dist().is_none() || prob_denovo.exp() < 0.95 {
            // no denovo variant, skip
            return None;
        }

        let mut vaf_dist = sample_info.vaf_dist().as_ref().unwrap().iter().map(|(vaf, density)| (*vaf, *density)).collect();

        Some(VariantObservation {
            vaf_dist,
            prob_denovo,
            max_posterior_vaf: sample_info.allelefreq_estimate()
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
                inf_density.ln_add_exp(LogProb(((sup_density.exp() - inf_density.exp()) / (**supremum - **infimum)).ln() + (vaf - *infimum).ln()))
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

#[derive(Serialize)]
struct EmpiricalVAFDistSpecEntry {
    vaf: f64,
    count: usize,
}

#[derive(Clone, Debug)]
struct EmpiricalVAFDist {
    inner: BTreeMap<AlleleFreq, usize>,
    total: usize,
    outliers: Vec<AlleleFreq>,
    bin_width: AlleleFreq,
}

impl EmpiricalVAFDist {
    fn new(variant_observations: &[VariantObservation]) -> Self {
        let mut inner = BTreeMap::default();
        let mut total = 0;
        for obs in variant_observations {
            let bin = AlleleFreq((*obs.max_posterior_vaf * 100.0).floor() / 100.0);
            *inner.entry(bin).or_insert(0) += 1;
            total += 1;
        }
        let dist = EmpiricalVAFDist {
            inner,
            total,
            outliers: Vec::new(),
            bin_width: AlleleFreq(0.1),
        };

        dist.strip_outliers()
    }

    fn as_json(&self) -> serde_json::Value {
        serde_json::Value::Array(self.inner.iter().map(|(vaf, count)| {
            json!({
                "vaf": *vaf,
                "count": *count
            })
        }).collect_vec())
    }

    fn bins<'a>(&'a self) -> impl Iterator<Item=f64> + 'a {
        self.inner.keys().map(|vaf| **vaf)
    }

    fn max_value(&self) -> AlleleFreq {
        self.inner.keys().last().cloned().unwrap_or(AlleleFreq(0.0))
    }

    fn strip_outliers(&self) -> Self {
        let mut stripped = self.clone();
        let mut outliers = stripped.outliers.clone();
        stripped.inner.retain(|vaf, count| {
            if *count < 2 {
                outliers.push(*vaf);
                false
            } else {
                *count -= 1;
                true
            }
        });

        stripped.outliers = outliers;

        stripped
    }

    fn scale_to(&self, max_value: AlleleFreq) -> Self {
        if *self.max_value() == 0.0 {
            // distribution is already empty
            return self.clone();
        }
        let factor = max_value / self.max_value();

        EmpiricalVAFDist {
            inner: self.inner.iter().map(|(vaf, count)| (*vaf * factor, *count)).collect(),
            outliers: self.outliers.iter().map(|vaf| *vaf * factor).collect(),
            total: self.total,
            bin_width: self.bin_width * factor,
        }
    }

    fn pdf(&self, vaf: AlleleFreq) -> LogProb {
        let denominator = self.total as f64;
        let count = self.inner.get(&vaf).expect("bug: unexpected vaf query");
        LogProb::from(Prob::from(*count as f64 / denominator))
    }
}

pub(crate) struct ContaminationEstimator {
    variant_observations: Vec<VariantObservation>,
    output: Option<PathBuf>,
    output_plot: PathBuf,
}

impl ContaminationEstimator {
    pub(crate) fn new(output: Option<PathBuf>) -> Self {
        ContaminationEstimator {
            variant_observations: Vec::new(),
            output,
            output_plot: "plot.vl.json".to_owned().into(),
        }
    }

    fn likelihood_single_obs(&self, variant_observation: &VariantObservation, contamination_expected_dist: &EmpiricalVAFDist) -> LogProb {
        LogProb::ln_trapezoidal_integrate_grid_exp(
            |_, vaf| {
                variant_observation.pdf(AlleleFreq(vaf)) + contamination_expected_dist.pdf(AlleleFreq(vaf))
            },
            &contamination_expected_dist.bins().collect_vec()
        )
    }

    fn likelihood(&self, contamination: AlleleFreq, empirical_vaf_dist: &EmpiricalVAFDist) -> LogProb {
        let contamination_expected_dist = empirical_vaf_dist.scale_to(AlleleFreq(1.0) - contamination);

        self.variant_observations.iter().map(|obs| self.likelihood_single_obs(obs, &contamination_expected_dist)).sum()
    }

    fn calc_posterior<W: Write>(&self, mut writer: csv::Writer<W>) -> Result<()> {
        let empirical_vaf_dist = EmpiricalVAFDist::new(&self.variant_observations);
        let contaminations = linspace(0.0, 1.0, 100).collect_vec();
        let likelihoods = contaminations.iter().map(|contamination| self.likelihood(AlleleFreq(*contamination), &empirical_vaf_dist)).collect_vec();

        // calculate marginal probability
        let marginal = LogProb::ln_trapezoidal_integrate_grid_exp(
            |i, _| likelihoods[i],
            &contaminations,
        );

        let mut spec = serde_json::from_str(include_str!("../../templates/plots/contamination_estimation.json"))?;
        if let serde_json::Value::Object(ref mut spec) = spec {
            spec["datasets"]["empirical_vaf_dist"] = empirical_vaf_dist.as_json();
            spec["datasets"]["posterior_density"] = serde_json::Value::Array(contaminations.iter().zip(likelihoods.iter()).map(|(contamination, likelihood)| {
                json!({
                    "purity": 1.0 - *contamination,
                    "density": (likelihood - marginal).exp()
                })
            }).collect_vec());

            let mut outfile = File::create(&self.output_plot)?;
            serde_json::to_writer_pretty(outfile, &spec)?;
        } else {
            unreachable!();
        }

        // sort in descending order
        let sorted_idx = (0..contaminations.len())
            .sorted_by_key(|i| NotNan::new(*likelihoods[*i]).unwrap())
            .collect_vec();
        let contaminations = sorted_idx
            .iter()
            .rev()
            .map(|i| contaminations[*i])
            .collect_vec();
        let likelihoods = sorted_idx
            .iter()
            .rev()
            .map(|i| likelihoods[*i])
            .collect_vec();

        // write into table
        writer.write_record(&["contamination", "posterior density"])?;
        for (contamination, likelihood) in contaminations.iter().zip(likelihoods.iter()) {
            let posterior = *likelihood - marginal;
            writer.write_record(&[
                format!("{}", contamination),
                format!("{}", *Prob::from(posterior)),
            ])?;
        }
        Ok(())
    }
}

impl CallProcessor for ContaminationEstimator {
    fn setup<Pr: bayesian::model::Prior, CF: CandidateFilter>(
        &mut self,
        caller: &Caller<Pr, Self, CF>,
    ) -> Result<()> {
        Ok(())
    }

    fn process_call(&mut self, call: Call, sample_names: &grammar::SampleInfo<String>) -> Result<()> {
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
        let contaminant_pileup = &work_item.pileups().as_ref().unwrap()[sample_names.iter().position(|s| *s == "contaminant").unwrap()];
        let sample_pileup = &work_item.pileups().as_ref().unwrap()[sample_names.iter().position(|s| *s == "sample").unwrap()];
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
        ContaminationEstimator::new(output),
        ContaminationCandidateFilter::new(),
    )
}

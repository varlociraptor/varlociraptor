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
use rgsl::randist::binomial::binomial_pdf;

use crate::{
    calling::variants::{
        calling::{call_generic, CallProcessor, Caller, CandidateFilter, WorkItem},
        Call,
    },
    grammar,
    utils::PathMap,
    variants::model::AlleleFreq,
    utils::NUMERICAL_EPSILON,
};

#[derive(Hash, Eq, PartialEq, Clone, Debug)]
struct Event {
    contamination: AlleleFreq,
    ploidy: u32,
}

#[derive(Clone, Debug)]
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

#[derive(new)]
struct Prior {
    histology_contamination: AlleleFreq,
    considered_cells: u64,
}

impl bio::stats::bayesian::model::Prior  for Prior {
    type Event = Event;

    fn compute(&self, event: &Self::Event) -> LogProb {
        // TODO: add pathology based prior
        LogProb::ln_one()
    }
}

#[derive(new)]
struct Likelihood {
    somatic_mutation_rate: f64,
    heterozygosity: f64,
    genome_size: f64,
}

impl bio::stats::bayesian::model::Likelihood for Likelihood {
    type Event = Event;
    type Data = Vec<VariantObservation>;

    fn compute(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        payload: &mut ()
    ) -> LogProb {
        let purity = AlleleFreq(1.0) - event.contamination;
        data.iter().map(|obs| {
            if *purity == 0.0 {
                // there cannot be any denovo somatic mutation
                return obs.prob_denovo.ln_one_minus_exp();
            }

            // clonal somatic mutations
            let prob_clonal = LogProb::ln_sum_exp(&(1..=event.ploidy).map(|m| {
                let prior = LogProb::from(Prob(self.heterozygosity / m as f64));
                let expected_vaf = AlleleFreq(m as f64 / event.ploidy as f64) * purity;
                obs.pdf(expected_vaf) + prior
            }).collect_vec());


            // subclonal somatic mutations
            let density = |_, vaf: f64| {
                let prior = LogProb((self.somatic_mutation_rate.ln()
                        - (2.0 * vaf.ln() + self.genome_size.ln())));
                let expected_vaf = AlleleFreq(vaf) * purity;
                obs.pdf(expected_vaf) + prior
            };
            let prob_subclonal = LogProb::ln_trapezoidal_integrate_grid_exp(
                density,
                &obs.vaf_dist.keys().filter_map(|vaf| if **vaf > 0.0 { Some(**vaf) } else { None }).collect_vec()
            );

            prob_clonal.ln_add_exp(prob_subclonal)
        }).sum()
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
        joint_prob: &mut F
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
        joint_prob: &mut F
    ) -> LogProb {
        LogProb::ln_sum_exp(
            &(1..8).map(|ploidy| {
                let density = |_, contamination| {
                    let event = Event { ploidy, contamination: AlleleFreq(contamination) };
                    joint_prob(&event, data)
                };
                LogProb::ln_simpsons_integrate_exp(
                    density,
                    0.0,
                    1.0,
                    51,
                )
            }).collect_vec()
        )
    }
}


#[derive(Clone, Debug)]
struct VAFDist {
    histogram: BTreeMap<AlleleFreq, usize>,
    total: usize,
}

impl VAFDist {
    fn new(variant_observations: &[VariantObservation]) -> Self {
        let mut histogram = BTreeMap::default();
        let mut total = 0;
        for obs in variant_observations {
            let bin = AlleleFreq((*obs.max_posterior_vaf * 100.0).floor() / 100.0);
            *histogram.entry(bin).or_insert(0) += 1;
            total += 1;
        }
        VAFDist {
            histogram,
            total,
        }
    }

    fn as_json(&self) -> serde_json::Value {
        serde_json::Value::Array(self.histogram.iter().map(|(vaf, count)| {
            json!({
                "vaf": *vaf,
                "count": *count
            })
        }).collect_vec())
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

    fn calc_posterior<W: Write>(&self, mut writer: csv::Writer<W>) -> Result<()> {
        let vaf_dist = VAFDist::new(&self.variant_observations);

        let prior = Prior::new(AlleleFreq(0.9), 100); // TODO parse from CLI
        let likelihood = Likelihood::new(1e-6, 0.001, 3.5e9);
        let model = bio::stats::bayesian::Model::new(likelihood, prior, Posterior::new());

        let model_instance = model.compute_from_marginal(&Marginal::new(), &self.variant_observations);

        let mut spec = serde_json::from_str(include_str!("../../templates/plots/contamination_estimation.json"))?;
        if let serde_json::Value::Object(ref mut spec) = spec {
            spec["datasets"]["empirical_vaf_dist"] = vaf_dist.as_json();
            spec["datasets"]["posterior_density"] = serde_json::Value::Array(model_instance.event_posteriors().map(|(event, density)| {
                json!({
                    "purity": 1.0 - *event.contamination,
                    "ploidy": event.ploidy,
                    "density": *Prob::from(density),
                })
            }).collect_vec());

            let mut outfile = File::create(&self.output_plot)?;
            serde_json::to_writer_pretty(outfile, &spec)?;
        } else {
            unreachable!();
        }

        // write into table
        writer.write_record(&["ploidy", "contamination", "posterior density"])?;
        for (event, density) in model_instance.event_posteriors() {
            writer.write_record(&[
                format!("{}", event.ploidy),
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

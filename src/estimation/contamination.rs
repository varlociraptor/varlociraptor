use std::{
    convert::TryFrom,
    io::{stdout, Write},
    path::PathBuf,
};

use anyhow::Result;
use bio::stats::{bayesian, LogProb, Prob};
use csv::WriterBuilder;
use itertools::Itertools;
use itertools_num::linspace;
use ordered_float::NotNan;

use crate::{
    calling::variants::{
        calling::{call_generic, CallProcessor, Caller, CandidateFilter, WorkItem},
        Call,
    },
    grammar,
    utils::PathMap,
    variants::model::AlleleFreq,
};

pub(crate) struct ContaminationEstimator {
    likelihoods: Vec<LogProb>,
    contaminations: Vec<f64>,
    output: Option<PathBuf>,
}

impl ContaminationEstimator {
    pub(crate) fn new(output: Option<PathBuf>) -> Self {
        let contaminations = linspace(0.0, 1.0, 100).collect();
        ContaminationEstimator {
            likelihoods: vec![LogProb::ln_one(); 100],
            contaminations,
            output,
        }
    }

    fn write_posterior<W: Write>(&self, mut writer: csv::Writer<W>) -> Result<()> {
        // calculate marginal probability
        let marginal = LogProb::ln_trapezoidal_integrate_grid_exp(
            |i, _| self.likelihoods[i],
            &self.contaminations,
        );

        // sort in descending order
        let sorted_idx = (0..self.contaminations.len())
            .sorted_by_key(|i| NotNan::new(*self.likelihoods[*i]).unwrap())
            .collect_vec();
        let contaminations = sorted_idx
            .iter()
            .rev()
            .map(|i| self.contaminations[*i])
            .collect_vec();
        let likelihoods = sorted_idx
            .iter()
            .rev()
            .map(|i| self.likelihoods[*i])
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
        let sample_info = call.variant().as_ref().unwrap().sample_info()[sample_names.iter().position(|s| *s == "sample").unwrap()]
            .as_ref()
            .unwrap();
        let prob_denovo = call.variant().as_ref().unwrap().event_probs().as_ref().unwrap().get("denovo").unwrap();
        let prob_not_denovo = prob_denovo.ln_one_minus_exp();

        if sample_info.vaf_dist().is_none() {
            // variant is most likely an artifact, skip
            return Ok(());
        }

        if prob_denovo.exp() >= 0.95 {
            //dbg!((call.chrom(), call.pos(), sample_info.allelefreq_estimate(), call.variant().as_ref().unwrap().event_probs()));
            println!("{}: {}", call.pos(), sample_info.allelefreq_estimate());
        } else {
            return Ok(());
        }

        for (contamination, factor) in self.contaminations.iter().zip(self.likelihoods.iter_mut()) {
            let mut vaf_dist = sample_info.vaf_dist().as_ref().unwrap().clone();
            // calculate likelihood factor for this variant:
            // the probability that the VAF is <=(1-contamination)=purity
            let purity = AlleleFreq(1.0 - *contamination);

            let mut vafs: Vec<_> = vaf_dist.keys().cloned().sorted().collect();
            let grid_points = match vafs.binary_search(&purity) {
                Ok(i) => vafs[..=i].iter().map(|vaf| **vaf).collect(),
                Err(i) => {
                    // i is the position where purity has to be inserted, i.e. the supremum
                    let mut grid_points: Vec<_> = vafs[..i].iter().map(|vaf| **vaf).collect();
                    grid_points.push(*purity);
                    if vafs.len() > i && i > 0 {
                        grid_points.push(*vafs[i]);
                        let infimum = vafs[i - 1];
                        let supremum = vafs[i];

                        // interpolate and insert into vaf_dist
                        let factor = (purity - infimum) / (purity - supremum);
                        let interpolated_prob = LogProb(
                            factor.ln()
                                + (vaf_dist.get(&supremum).unwrap().exp()
                                    - vaf_dist.get(&infimum).unwrap().exp())
                                .ln(),
                        );
                        vaf_dist.insert(purity, interpolated_prob);
                    } else if i > 0 {
                        // nothing to interpolate, just repeat the infimum prob
                        vaf_dist.insert(purity, *vaf_dist.get(&vafs[i - 1]).unwrap());
                    } else {
                        // purity comes at the beginning, repeat supremum prob
                        vaf_dist.insert(purity, *vaf_dist.get(&vafs[0]).unwrap());
                    }

                    grid_points
                }
            };

            *factor += (LogProb::ln_trapezoidal_integrate_grid_exp(
                |_, vaf| vaf_dist.get(&AlleleFreq(vaf)).cloned().unwrap(),
                &grid_points,
            ) + prob_denovo).ln_add_exp(prob_not_denovo);
        }

        Ok(())
    }

    fn finalize(&mut self) -> Result<()> {
        if let Some(ref path) = self.output {
            self.write_posterior(WriterBuilder::new().delimiter(b'\t').from_path(path)?)
        } else {
            self.write_posterior(WriterBuilder::new().delimiter(b'\t').from_writer(stdout()))
        }
    }
}

#[derive(new)]
struct ContaminationCandidateFilter;

impl CandidateFilter for ContaminationCandidateFilter {
    fn filter(&self, work_item: &WorkItem, sample_names: &grammar::SampleInfo<String>) -> bool {
        // only keep variants where all reads of the contaminant are REF and that are SNVs
        let contaminant_pileup = &work_item.pileups().as_ref().unwrap()[sample_names.iter().position(|s| *s == "contaminant").unwrap()];
        work_item.snv().is_some()
            && (contaminant_pileup.read_observations().len() <= 10
            || contaminant_pileup
                .read_observations()
                .iter()
                .all(|obs| obs.is_ref_support()))
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

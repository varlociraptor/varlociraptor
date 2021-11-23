use anyhow::Result;
use bio::stats::probs::{LogProb, Prob};
use derive_builder::Builder;
use hdf5;
use ordered_float::NotNaN;
use rust_htslib::bcf;
use serde::Serialize;
use std::collections::HashMap;
use std::io;
use std::path::PathBuf;

use bio::stats::bayesian::model::Model;

use crate::haplotypes::model::{Data, HaplotypeFractions, Likelihood, Marginal, Posterior, Prior};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    hdf5_reader: hdf5::File,
    vcf_reader: bcf::Reader,
    min_norm_counts: f64,
    outcsv: Option<PathBuf>,
}

#[derive(Serialize)]
struct Output {
    haplotype_a: NotNaN<f64>,
    haplotype_b: NotNaN<f64>,
    haplotype_c: NotNaN<f64>,
    posterior: f64,
    odds_ratio: f64,
}

impl Caller {
    pub fn call(&self) -> Result<()> {
        // Step 1: obtain kallisto estimates.
        let kallisto_estimates = KallistoEstimates::new(&self.hdf5_reader, self.min_norm_counts)?;

        // Step 2: setup model.
        let model = Model::new(Likelihood::new(), Prior::new(), Posterior::new());

        //let universe = HaplotypeFractions::likely(&kallisto_estimates);
        let data = Data::new(kallisto_estimates.values().cloned().collect());

        // Step 3: calculate posteriors.
        //let m = model.compute(universe, &data);
        let n_of_haplotypes = 3;
        let m = model.compute_from_marginal(&Marginal::new(n_of_haplotypes), &data);

        // Step 4: print TSV table with results
        // TODO use csv crate
        // Columns: posterior_prob, haplotype_a, haplotype_b, haplotype_c, ...
        // with each column after the first showing the fraction of the respective haplotype
        let mut posterior = m.event_posteriors();
        let mut wtr = csv::Writer::from_path(self.outcsv.as_ref().unwrap())?;
        let mut haplotype_fractions = Vec::new();
        let best = posterior.next().unwrap();
        let best_density = best.1 .0.exp();

        for haplotype in 0..n_of_haplotypes {
            let frequencies = best.0;
            haplotype_fractions.push(frequencies[haplotype]);
        }
        wtr.serialize(Output {
            haplotype_a: haplotype_fractions[0],
            haplotype_b: haplotype_fractions[1],
            haplotype_c: haplotype_fractions[2],
            posterior: best_density,
            odds_ratio: 1.0, //for the best record
        })?;

        for (haplotype_frequencies, logprob) in posterior {
            let density = logprob.0.exp();
            let odds = density / best_density;
            let mut haplotype_fractions = Vec::new();
            for frequency in haplotype_frequencies.iter() {
                haplotype_fractions.push(frequency);
            }
            wtr.serialize(Output {
                haplotype_a: *haplotype_fractions[0],
                haplotype_b: *haplotype_fractions[1],
                haplotype_c: *haplotype_fractions[2],
                posterior: density,
                odds_ratio: odds,
            })?;
        }
        Ok(())
    }
}

#[derive(Derefable, Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct Haplotype(#[deref] String);

#[derive(Debug, Clone)]
pub(crate) struct KallistoEstimate {
    pub count: f64,
    pub dispersion: f64,
}

#[derive(Debug, Clone, Derefable)]
pub(crate) struct KallistoEstimates(#[deref] HashMap<Haplotype, KallistoEstimate>);

impl KallistoEstimates {
    /// Generate new instance.
    pub(crate) fn new(hdf5_reader: &hdf5::File, min_norm_counts: f64) -> Result<Self> {
        let seqnames = Self::filter_seqnames(hdf5_reader, min_norm_counts)?;

        let ids = hdf5_reader
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::FixedAscii<255>>()?;
        let num_bootstraps = hdf5_reader.dataset("aux/num_bootstrap")?.read_1d::<i32>()?;
        let seq_length = hdf5_reader.dataset("aux/lengths")?.read_1d::<f64>()?;

        let mut estimates = HashMap::new();

        for seqname in seqnames {
            let index = ids.iter().position(|x| x.as_str() == seqname).unwrap();
            let mut bootstraps = Vec::new();
            for i in 0..num_bootstraps[0] {
                let dataset = hdf5_reader.dataset(&format!("bootstrap/bs{i}", i = i))?;
                let est_counts = dataset.read_1d::<f64>()?;
                let norm_counts = est_counts / &seq_length;
                let norm_counts = norm_counts[index];
                bootstraps.push(norm_counts);
            }
            //mean
            let sum = bootstraps.iter().sum::<f64>();
            let count = bootstraps.len();
            let m = sum / count as f64;

            //std dev
            let variance = bootstraps
                .iter()
                .map(|value| {
                    let diff = m - (*value as f64);
                    diff * diff
                })
                .sum::<f64>()
                / count as f64;
            let std = variance.sqrt();
            let t = std / m;

            //retrieval of mle
            let mle_dataset = hdf5_reader.dataset("est_counts")?.read_1d::<f64>()?;
            let mle_norm = mle_dataset / &seq_length; //normalized mle counts by length
            let m = mle_norm[index];

            estimates.insert(
                Haplotype(seqname.clone()),
                KallistoEstimate {
                    dispersion: t,
                    count: m,
                },
            );
        }

        Ok(KallistoEstimates(estimates))
    }

    //Return a vector of filtered seqnames according to --min-norm-counts.
    fn filter_seqnames(hdf5_reader: &hdf5::File, min_norm_counts: f64) -> Result<Vec<String>> {
        let est_counts = hdf5_reader.dataset("est_counts")?.read_1d::<f64>()?;
        let seq_length = hdf5_reader.dataset("aux/lengths")?.read_1d::<f64>()?; //these two variables arrays have the same length.
        let norm_counts = est_counts / seq_length;
        let mut indices = Vec::new();
        for (i, num) in norm_counts.iter().enumerate() {
            if num > &min_norm_counts {
                indices.push(i);
            }
        }
        let ids = hdf5_reader
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::FixedAscii<255>>()?;

        let mut filtered: Vec<String> = Vec::new();
        for i in indices {
            filtered.push(ids[i].to_string());
        }
        Ok(filtered)
    }
}

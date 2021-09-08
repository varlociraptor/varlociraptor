use anyhow::Result;
use bio::stats::probs::LogProb;
use derive_builder::Builder;
use hdf5;
use kernel_density;
use ordered_float::OrderedFloat;
use rust_htslib::bcf;
use serde_json::json;
use statrs::function::beta::ln_beta;
use std::fs::File;
use std::mem;
use std::path::PathBuf;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct Caller {
    hdf5_reader: hdf5::File,
    vcf_reader: bcf::Reader,
    min_norm_counts: f64,
}

#[derive(Debug, Deserialize, Serialize)]
pub(crate) struct ECDF {
    seqname: String,
    bootstraps: Vec<f64>,
    ecdf: Vec<f64>,
}

#[derive(Debug, Deserialize, Serialize)]
pub(crate) struct PlotEcdf {
    bootstrap: f64,
    ecdf: f64,
}

impl Caller {
    //Below function outputs a vector of filtered seqnames according to --min-norm-counts.
    pub(crate) fn filter_seqnames(&self) -> Result<Vec<String>> {
        let hdf5 = &self.hdf5_reader;
        let est_counts = hdf5.dataset("est_counts")?.read_1d::<f64>()?;
        let seq_length = hdf5.dataset("aux/lengths")?.read_1d::<f64>()?; //these two variables arrays have the same length.
        let norm_counts = est_counts / seq_length;
        let mut indices = Vec::new();
        for (i, num) in norm_counts.iter().enumerate() {
            if num > &self.min_norm_counts {
                indices.push(i);
            }
        }
        //let len = hdf5.dataset("aux/ids")?.dtype().unwrap().size();
        let ids = hdf5
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::FixedAscii<255>>()?;

        let mut filtered: Vec<String> = Vec::new();
        for i in indices {
            filtered.push(ids[i].to_string());
        }
        Ok(filtered)
    }

    //Below function outputs an ECDF struct for one haplotype.
    pub(crate) fn cdf(&self, seqname: String) -> Result<ECDF> {
        let hdf5 = &self.hdf5_reader;
        let ids = hdf5
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::FixedAscii<255>>()?;
        let index = ids.iter().position(|x| x.as_str() == seqname).unwrap();
        let num_bootstraps = hdf5.dataset("aux/num_bootstrap")?.read_1d::<i32>()?;
        let seq_length = hdf5.dataset("aux/lengths")?.read_1d::<f64>()?;
        let mut bootstraps = Vec::new();
        for i in 0..num_bootstraps[0] {
            let dataset = hdf5.dataset(&format!("bootstrap/bs{i}", i = i))?;
            let est_counts = dataset.read_1d::<f64>()?;
            let norm_counts = est_counts / &seq_length;
            let norm_counts = norm_counts[index];
            bootstraps.push(norm_counts);
        }
        let bootstraps_ordered_float: Vec<OrderedFloat<f64>> =
            bootstraps.iter().map(|&x| OrderedFloat(x)).collect();
        let ecdf = kernel_density::ecdf::Ecdf::new(&bootstraps_ordered_float); //the trait Ord is not implemented for f64!
        let ecdf: Vec<f64> = bootstraps_ordered_float
            .iter()
            .map(|&x| ecdf.value(x))
            .collect();
        Ok(ECDF {
            seqname,
            bootstraps,
            ecdf,
        })
    }

    pub(crate) fn app_neg_binom(&self, seqname: String) -> Result<()> {
        let hdf5 = &self.hdf5_reader;
        let ids = hdf5
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::FixedAscii<255>>()?;
        let index = ids.iter().position(|x| x.as_str() == seqname).unwrap();
        let num_bootstraps = hdf5.dataset("aux/num_bootstrap")?.read_1d::<i32>()?;
        let seq_length = hdf5.dataset("aux/lengths")?.read_1d::<f64>()?;
        let mut bootstraps = Vec::new();
        for i in 0..num_bootstraps[0] {
            let dataset = hdf5.dataset(&format!("bootstrap/bs{i}", i = i))?;
            let est_counts = dataset.read_1d::<f64>()?;
            let norm_counts = est_counts / &seq_length;
            let norm_counts = norm_counts[index];
            bootstraps.push(norm_counts);
        }

        //mean
        let sum = bootstraps.iter().sum::<f64>();
        let count = bootstraps.len();
        let mean = sum / count as f64;

        //std dev
        let variance = bootstraps
            .iter()
            .map(|value| {
                let diff = mean - (*value as f64);
                diff * diff
            })
            .sum::<f64>()
            / count as f64;
        let std = variance.sqrt();

        let theta = std / mean;

        let mle_dataset = hdf5.dataset("est_counts")?.read_1d::<f64>()?;
        let mle_norm = mle_dataset / &seq_length; //normalized mle counts by length
        let mle = mle_norm[index];

        let neg_bin = neg_binom(mle, mean, theta);
        Ok(())
    }

    pub(crate) fn call(&self) {
        todo!()
    }
}

impl ECDF {
    pub(crate) fn plot_qc(&self, mut qc_plot: PathBuf) -> Result<()> {
        let json = include_str!("../../templates/plots/qc_cdf.json");
        let mut blueprint: serde_json::Value = serde_json::from_str(json)?;

        let mut plot_data = Vec::new();

        for num in 0..self.bootstraps.to_vec().len() {
            plot_data.push(PlotEcdf {
                bootstrap: self.bootstraps[num],
                ecdf: self.ecdf[num],
            })
        }
        let plot_data = json!(plot_data);
        blueprint["data"]["values"] = plot_data;

        qc_plot.push(self.seqname.clone() + &".json".to_string());
        let file = File::create(qc_plot)?;
        serde_json::to_writer(file, &blueprint);
        Ok(())
    }
}

fn neg_binom(x: f64, mu: f64, theta: f64) -> LogProb {
    let n = 1.0 / theta;
    let p = n / (n + mu);
    let mut p1 = if n > 0.0 { n * p.ln() } else { 0.0 };
    let mut p2 = if x > 0.0 { x * (1.0 - p).ln() } else { 0.0 };
    let b = ln_beta(x + 1.0, n);

    if p1 < p2 {
        mem::swap(&mut p1, &mut p2);
    }
    LogProb((p1 - b + p2) - (x + n).ln())
}

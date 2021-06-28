use anyhow::Result;
use derive_builder::Builder;
use hdf5;
use kernel_density;
use rust_htslib::bcf;
use serde_json::json;
use std::fs::File;
use std::path::PathBuf;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct Caller {
    hdf5_reader: hdf5::File,
    vcf_reader: bcf::Reader,
    min_norm_counts: u64,
    qc_plot: Option<PathBuf>,
}

pub(crate) struct ECDF {
    seqname: String,
    bootstraps: Vec<u64>,
    ecdf: Vec<f64>,
}

impl Caller {
    //Below function outputs a vector of filtered seqnames according to --min-norm-counts.
    pub(crate) fn filter_seqnames(&self) -> Result<Vec<String>> {
        let hdf5 = &self.hdf5_reader;
        let est_counts = hdf5.dataset("/est_counts")?.read_1d::<u64>()?;
        let seq_length = hdf5.dataset("/aux/lengths")?.read_1d::<u64>()?; //these two variables arrays have the same length.

        let norm_counts = est_counts / seq_length;

        let mut indices = Vec::new();
        for (i, num) in norm_counts.iter().enumerate() {
            if num > &self.min_norm_counts {
                indices.push(i);
            }
        }
        let ids = hdf5
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::VarLenAscii>()?;

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
            .read_1d::<hdf5::types::VarLenAscii>()?;
        let index = ids.iter().position(|x| x.as_str() == seqname).unwrap();
        let num_bootstraps = hdf5.dataset("/aux/num_bootstrap")?.read_scalar()?;
        let seq_length = hdf5.dataset("/aux/lengths")?.read_1d::<u64>()?;
        let mut bootstraps = Vec::new();
        for i in 0..num_bootstraps {
            let dataset = hdf5.dataset(&format!("bootstrap/bs{i}", i = i))?;
            let est_counts = dataset.read_1d::<u64>()?;
            let norm_counts = est_counts / &seq_length;
            let norm_counts = norm_counts[index];
            bootstraps.push(norm_counts);
        }
        let ecdf = kernel_density::ecdf::Ecdf::new(&bootstraps);
        let ecdf: Vec<f64> = bootstraps.iter().map(|&x| ecdf.value(x)).collect();
        Ok(ECDF {
            seqname,
            bootstraps,
            ecdf,
        })
    }

    pub(crate) fn call(&self) {
        todo!()
    }
}

impl ECDF {
    pub(crate) fn plot_qc(&self, qc_plot: Option<PathBuf>) -> Result<()> {
        let blueprint = "../../templates/plots/qc_cdf.json";
        let mut blueprint: serde_json::Value = serde_json::from_str(blueprint)?;

        let bootstraps = &self.bootstraps;
        let ecdf = &self.ecdf;

        let data = json!([bootstraps, ecdf]);
        blueprint["data"]["values"] = data;

        serde_json::to_writer(&File::create(qc_plot.as_ref().unwrap())?, &blueprint)?;
        Ok(())
    }
}

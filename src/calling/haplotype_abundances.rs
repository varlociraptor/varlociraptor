use anyhow::Result;
use derive_builder::Builder;
use hdf5;
use kernel_density;
use rust_htslib::bcf;


#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct Caller {
    hdf5_reader: hdf5::File,
    vcf_reader: bcf::Reader,
}

impl Caller {
    pub(crate) fn cdf(&self) -> Result<()> {
        let hdf5 = &self.hdf5_reader;
        let group = hdf5.group("bootstrap")?;
        let num_bootstraps = hdf5.dataset("/aux/num_bootstrap")?.read_scalar()?;
        let mut bootstraps = Vec::new();
        for i in 0..num_bootstraps {
            let dataset = group.dataset(&format!("bs{i}", i = i))?;
            let tpm = dataset.read_1d::<u64>()?;
            bootstraps.extend(tpm.to_vec());
        }
        let n_values = hdf5.dataset("/aux/lengths")?.read_scalar()?;
        let array = ndarray::Array2::from_shape_vec((num_bootstraps, n_values), bootstraps)?;
        let array = array.t();
        let cdf: Vec<_> = array
            .outer_iter()
            .map(|x| kernel_density::ecdf::Ecdf::new(&x.to_vec()))
            .collect();
        Ok(())
    }
    pub(crate) fn call(&self) {
        todo!()
    }
}

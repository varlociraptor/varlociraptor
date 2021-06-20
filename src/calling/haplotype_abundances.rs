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

pub(crate) struct ECDF {
    seqname: String,
    bootstraps: Vec<u64>,
    ecdf: kernel_density::ecdf::Ecdf<u64>,
}

impl Caller {
    pub(crate) fn cdf(&self, seqname: String) -> Result<ECDF> {
        let hdf5 = &self.hdf5_reader;
        let ids = hdf5
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::VarLenAscii>()?;
        let index = ids.iter().position(|x| x.as_str() == seqname).unwrap();
        let num_bootstraps = hdf5.dataset("/aux/num_bootstrap")?.read_scalar()?;
        let mut bootstraps = Vec::new();
        for i in 0..num_bootstraps {
            let dataset = hdf5.dataset(&format!("bootstrap/bs{i}", i = i))?;
            let tpm = dataset.read_1d::<u64>()?;
            let tst = tpm[index];
            bootstraps.push(tst);
        }
        let ecdf = kernel_density::ecdf::Ecdf::new(&bootstraps);
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

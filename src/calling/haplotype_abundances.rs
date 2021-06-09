use hdf5;
use rust_htslib::bcf;
use derive_builder::Builder;
use anyhow::Result;
use std::path::Path;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct Caller<'a> {
    hdf5_reader: hdf5::Reader<'a>,
    vcf_reader: bcf::Reader,
}

impl<'a> Caller<'a> {
    pub(crate) fn call(&self) {
        todo!()
    }
}


// impl<'a> CallerBuilder<'a> {
//     pub(crate) fn new(hdf5_reader: hdf5::Reader, vcf_reader: bcf::Reader) -> CallerBuilder<'a> {
//         CallerBuilder {
//             hdf5_reader: hdf5::File::open(hdf5_reader.as_ref()),
//             vcf_reader: bcf::Reader::from_path(vcf_reader)
//         }
//     }
// }

// impl<'a> Caller<'a> {
//     pub(crate) fn new<P: AsRef<Path>>(hdf5_path: P, vcf_path: P) -> Result<Self> {
//         let hdf5_reader = hdf5::File::open(hdf5_path.as_ref())?;
//         let vcf_reader = bcf::Reader::from_path(vcf_path)?;
        
//         Ok (Self {
//             hdf5_reader,
//             vcf_reader,
//         })
//     }
//     pub(crate) fn call(&mut self) {
//         todo!()
//     }
//   }

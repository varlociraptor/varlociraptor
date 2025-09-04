use std::{collections::HashMap, path::Path};

use anyhow::Result;
use rust_htslib::bcf::{self, Read};

use crate::utils;
use crate::variants::model::HaplotypeIdentifier;

#[derive(Default, Debug)]
pub(crate) struct HaplotypeFeatureIndex {
    last_records: HashMap<HaplotypeIdentifier, usize>,
}

impl HaplotypeFeatureIndex {
    pub(crate) fn new<P: AsRef<Path>>(inbcf: P) -> Result<Self> {
        let mut bcf_reader = bcf::Reader::from_path(inbcf)?;
        if !utils::is_haplotype_bcf(&bcf_reader) {
            return Ok(HaplotypeFeatureIndex::default());
        }

        let mut last_records = HashMap::new();

        let mut i = 0;
        loop {
            let mut record = bcf_reader.empty_record();
            match bcf_reader.read(&mut record) {
                None => return Ok(HaplotypeFeatureIndex { last_records }),
                Some(res) => res?,
            }

            if let Some(identifier) = HaplotypeIdentifier::from(&mut record)? {
                last_records.insert(identifier, i);
            }

            i += 1;
        }
    }

    pub(crate) fn last_record_index(
        &self,
        haplotype_identifier: &HaplotypeIdentifier,
    ) -> Option<usize> {
        self.last_records.get(haplotype_identifier).cloned()
    }
}

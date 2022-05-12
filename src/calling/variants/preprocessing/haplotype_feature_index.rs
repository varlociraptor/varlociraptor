use std::{collections::HashMap, convert::TryFrom, path::Path};

use anyhow::Result;
use rust_htslib::bcf::{
    self,
    record::Numeric,
    record::{Genotype, GenotypeAllele},
    Read,
};

use crate::utils;

#[derive(Hash, Debug, Eq, PartialEq)]
pub(crate) enum HaplotypeIdentifier {
    Event(Vec<u8>),
    Phase { phase_set: u32, genotype: Genotype },
}

impl HaplotypeIdentifier {
    pub(crate) fn from(record: &mut bcf::Record) -> Result<Option<Self>> {
        if utils::is_bnd(record)? {
            // TODO support records without EVENT tag.
            if let Ok(Some(event)) = record.info(b"EVENT").string() {
                let event = event[0];
                return Ok(Some(HaplotypeIdentifier::Event(event.to_owned())));
            }
        } else {
            if let Ok(ps_values) = record.format(b"PS").integer() {
                let phase_set = ps_values[0][0];
                if !phase_set.is_missing() {
                    // phase set definition found
                    if let Ok(genotypes) = record.genotypes() {
                        let genotype = genotypes.get(0);
                        if genotype.iter().any(|allele| match allele {
                            GenotypeAllele::Phased(allele_idx) if *allele_idx > 0 => true,
                            _ => false,
                        }) {
                            return Ok(Some(HaplotypeIdentifier::Phase {
                                phase_set: phase_set as u32,
                                genotype,
                            }));
                        }
                    }
                }
            }
        }
        Ok(None)
    }
}

#[derive(Default, Debug)]
pub(crate) struct HaplotypeFeatureIndex {
    last_records: HashMap<HaplotypeIdentifier, usize>,
}

impl HaplotypeFeatureIndex {
    pub(crate) fn new<P: AsRef<Path>>(inbcf: P) -> Result<Self> {
        let mut bcf_reader = bcf::Reader::from_path(inbcf)?;
        if !utils::is_sv_bcf(&bcf_reader) {
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

use std::collections::{HashMap, HashSet};

use anyhow::Result;
use linear_map::LinearMap;
use rust_htslib::bcf::Read;
use rust_htslib::bcf::{self, header::TagType};

#[derive(Debug, Clone, Default)]
pub(crate) struct AuxInfo {
    values: HashMap<Vec<u8>, AuxInfoValue>,
}

impl AuxInfo {
    pub(crate) fn write(&self, record: &mut bcf::Record, omit: &HashSet<Vec<u8>>) -> Result<()> {
        for (field, value) in &self.values {
            if omit.contains(field) {
                continue;
            }
            match value {
                AuxInfoValue::Integer(values) => {
                    record.push_info_integer(field, values)?;
                }
                AuxInfoValue::Float(values) => {
                    record.push_info_float(field, values)?;
                }
                AuxInfoValue::String(values) => {
                    let values: Vec<&[u8]> = values.iter().map(|value| value.as_slice()).collect();
                    record.push_info_string(field, &values)?;
                }
                AuxInfoValue::Flag(value) => {
                    if *value {
                        record.push_info_flag(field)?;
                    }
                }
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub(crate) enum AuxInfoValue {
    Integer(Vec<i32>),
    Float(Vec<f32>),
    String(Vec<Vec<u8>>),
    Flag(bool),
}

impl AuxInfoValue {}

#[derive(Debug, Clone)]
pub(crate) struct AuxInfoCollector {
    fields: HashMap<Vec<u8>, TagType>,
    records: Vec<LinearMap<String, String>>,
}

impl AuxInfoCollector {
    pub(crate) fn new(fields: &[Vec<u8>], bcf_reader: &bcf::Reader) -> Result<Self> {
        let header = bcf_reader.header();

        let fields = fields
            .iter()
            .map(|field| Ok((field.to_owned(), header.info_type(field)?.0)))
            .collect::<Result<HashMap<Vec<u8>, TagType>>>()?;

        let records = header
            .header_records()
            .into_iter()
            .filter_map(|record| {
                if let bcf::header::HeaderRecord::Info { values, .. } = record {
                    if fields.contains_key(values.get("ID").unwrap().as_bytes()) {
                        Some(values)
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect();

        Ok(AuxInfoCollector { fields, records })
    }

    pub(crate) fn write_header_info(&self, header: &mut bcf::Header) {
        for record in &self.records {
            header.push_record(
                format!(
                    "##INFO=<ID={id},Number={number},Type={_type},Description={desc}>",
                    id = record.get("ID").unwrap(),
                    number = record.get("Number").unwrap(),
                    desc = record.get("Description").unwrap(),
                    _type = record.get("Type").unwrap(),
                )
                .as_bytes(),
            );
        }
    }

    pub(crate) fn collect(&self, record: &bcf::Record) -> Result<AuxInfo> {
        let mut aux_info = HashMap::new();
        for (field, tag_type) in &self.fields {
            match tag_type {
                TagType::Integer => {
                    let value = record.info(field).integer()?;
                    if let Some(value) = value {
                        aux_info.insert(field.to_owned(), AuxInfoValue::Integer(value.to_owned()));
                    }
                }
                TagType::Float => {
                    let value = record.info(field).float()?;
                    if let Some(value) = value {
                        aux_info.insert(field.to_owned(), AuxInfoValue::Float(value.to_owned()));
                    }
                }
                TagType::String => {
                    let value = record.info(field).string()?;
                    if let Some(value) = value {
                        aux_info.insert(
                            field.to_owned(),
                            AuxInfoValue::String(
                                value.iter().map(|value| value.to_vec()).collect(),
                            ),
                        );
                    }
                }
                TagType::Flag => {
                    let value = record.info(field).flag()?;
                    aux_info.insert(field.to_owned(), AuxInfoValue::Flag(value));
                }
            }
        }
        Ok(AuxInfo { values: aux_info })
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use rust_htslib::bcf;
    use std::path::Path;

    /// Create an `AuxInfoCollector` from an existing VCF file path.
    ///
    /// The VCF must already declare any fields in `fields` in its header.
    /// Pass empty `fields` slice for a no-op collector.
    ///
    /// # Arguments
    /// * `vcf_path` - Path to VCF file whose header declares the fields
    /// * `fields` - INFO field names to collect, empty for no propagation
    pub(crate) fn make_aux_collector(vcf_path: &Path, fields: &[&str]) -> AuxInfoCollector {
        let reader = bcf::Reader::from_path(vcf_path).unwrap();
        let fields_bytes: Vec<Vec<u8>> = fields.iter().map(|s| s.as_bytes().to_vec()).collect();
        AuxInfoCollector::new(&fields_bytes, &reader).unwrap()
    }
}

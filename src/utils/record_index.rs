use std::cell::RefCell;
use std::cmp;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::path::Path;
use std::rc::Rc;
use std::str;
use std::sync::Arc;

use anyhow::Result;
use bio::alphabets::dna;
use bio::stats::pairhmm::EmissionParameters;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};
use regex::Regex;
use rust_htslib::bam;
use rust_htslib::bcf::{self, Read};
use vec_map::VecMap;

use crate::errors::Error;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::model;
use crate::variants::sampling_bias::{ReadSamplingBias, SamplingBias};
use crate::variants::types::{
    AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocus, SingleLocusBuilder, Variant,
};
use crate::{default_emission, default_ref_base_emission};

#[derive(Default, Debug)]
pub(crate) struct RecordIndex {
    last_records: HashMap<Vec<u8>, usize>,
}

impl RecordIndex {
    pub(crate) fn new<P: AsRef<Path>>(inbcf: P) -> Result<Self> {
        let mut bcf_reader = bcf::Reader::from_path(inbcf)?;
        if !utils::is_sv_bcf(&bcf_reader) {
            return Ok(RecordIndex::default());
        }

        let mut last_records = HashMap::new();

        let mut i = 0;
        loop {
            let mut record = bcf_reader.empty_record();
            if !bcf_reader.read(&mut record)? {
                return Ok(RecordIndex { last_records });
            }

            if utils::is_bnd(&mut record)? {
                // TODO support records without EVENT tag.
                if let Ok(Some(event)) = record.info(b"EVENT").string() {
                    let event = event[0];
                    last_records.insert(event.to_owned(), i);
                }
            }

            i += 1;
        }
    }

    pub(crate) fn last_record_index(&self, event: &[u8]) -> Option<usize> {
        self.last_records.get(event).cloned()
    }
}
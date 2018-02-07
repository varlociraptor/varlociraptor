use std::str;
use std::collections::vec_deque;
use std::cmp;

use vec_map::VecMap;
use serde::Serialize;
use serde::ser::{Serializer, SerializeStruct};

use bio::stats::LogProb;
use rust_htslib::bam;
use rust_htslib::bam::record::CigarString;

use model::Variant;
use model::evidence;
use model::sample::RecordBuffer;


/// An observation for or against a variant.
#[derive(Clone)]
pub struct Observation {
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    pub prob_mapping: LogProb,
    /// Probability that the read/read-pair comes from the alternative allele.
    pub prob_alt: LogProb,
    /// Probability that the read/read-pair comes from the reference allele.
    pub prob_ref: LogProb,
    /// Probability of the read/read-pair given that it has been mismapped.
    pub prob_mismapped: LogProb,
    /// Read lengths
    pub left_read_len: u32,
    pub right_read_len: Option<u32>,
    /// Type of evidence.
    pub evidence: Evidence
}


impl Observation {
    pub fn is_alignment_evidence(&self) -> bool {
        if let Evidence::Alignment(_) = self.evidence {
            true
        } else {
            false
        }
    }
}


impl<'a> Serialize for Observation<'a> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: Serializer {
        let mut s = serializer.serialize_struct("Observation", 3)?;
        s.serialize_field("prob_mapping", &self.prob_mapping)?;
        s.serialize_field("prob_alt", &self.prob_alt)?;
        s.serialize_field("prob_ref", &self.prob_ref)?;
        s.serialize_field("evidence", &self.evidence)?;
        s.end()
    }
}


/// Types of evidence that lead to an observation.
/// The contained information is intended for debugging and will be printed together with
/// observations.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Evidence {
    /// Insert size of fragment
    InsertSize(String),
    /// Alignment of a single read
    Alignment(String)
}


impl Evidence {
    /// Create a dummy alignment.
    pub fn dummy_alignment() -> Self {
        Evidence::Alignment("Dummy-Alignment".to_owned())
    }

    /// Create dummy insert size evidence.
    pub fn dummy_insert_size(insert_size: u32) -> Self {
        Evidence::InsertSize(format!("insert-size={}", insert_size))
    }

    /// Create insert size evidence.
    pub fn insert_size(
        insert_size: u32,
        left: &CigarString,
        right: &CigarString,
        left_record: &bam::Record,
        right_record: &bam::Record,
        p_left_ref: LogProb,
        p_left_alt: LogProb,
        p_right_ref: LogProb,
        p_right_alt: LogProb,
        p_isize_ref: LogProb,
        p_isize_alt: LogProb
    ) -> Self {
        Evidence::InsertSize(format!(
            "left: cigar={} ({:e} vs {:e}), right: cigar={} ({:e} vs {:e}), insert-size={} ({:e} vs {:e}), qname={}, left: AS={:?}, XS={:?}, right: AS={:?}, XS={:?}",
            left, p_left_ref.exp(), p_left_alt.exp(),
            right, p_right_ref.exp(), p_right_alt.exp(),
            insert_size, p_isize_ref.exp(), p_isize_alt.exp(),
            str::from_utf8(left_record.qname()).unwrap(),
            left_record.aux(b"AS").map(|a| a.integer()),
            left_record.aux(b"XS").map(|a| a.integer()),
            right_record.aux(b"AS").map(|a| a.integer()),
            right_record.aux(b"XS").map(|a| a.integer())
        ))
    }

    /// Create alignment evidence.
    pub fn alignment(cigar: &CigarString, record: &bam::Record) -> Self {
        Evidence::Alignment(format!(
            "cigar={}, qname={}, AS={:?}, XS={:?}",
            cigar, str::from_utf8(record.qname()).unwrap(),
            record.aux(b"AS").map(|a| a.integer()),
            record.aux(b"XS").map(|a| a.integer())
        ))
    }
}


/// Data that is shared among all observations over a locus.
#[derive(Clone)]
pub struct Common {
    pub softclip_obs: Option<VecMap<u32>>,
    /// Average number of reads starting at any position in the region.
    pub coverage: f64,
    pub max_read_len: u32,
    pub enclosing_possible: bool
}


impl Common {
    pub fn new(records: &RecordBuffer, variant: &Variant) -> Self {
        let valid_records = || records.iter().filter(|rec| !rec.is_supplementary());

        // obtain maximum read len
        let max_read_len = valid_records().map(
            |rec| rec.seq().len()
        ).max().unwrap_or(0) as u32;

        // average number of reads starting at any position in the current window
        let coverage = valid_records().count() as f64 /
                       (records.window as f64 * 2.0);

        // determine if variant can be enclosed
        let enclosing_possible = if variant.is_snv() {
            true
        } else {
            let max_indel_cigar = valid_records().map(|rec| {
                evidence::max_indel(&rec.cigar())
            }).max().unwrap_or(0);
            variant.len() <= max_indel_cigar
        };

        let softclip_obs = if variant.is_snv() {
            None
        } else {
            let mut obs = VecMap::new();
            for rec in valid_records() {
                let cigar = rec.cigar();
                let s = cmp::max(
                    evidence::Clips::trailing(&cigar).soft(),
                    evidence::Clips::leading(&cigar).soft()
                );
                obs.entry(&s).or_insert(0) += 1;
            }

            Some(obs)
        };


        Common {
            softclip_obs: softclip_obs,
            coverage: coverage,
            max_read_len: max_read_len,
            enclosing_possible: enclosing_possible
        }
    }

    pub fn collect_softclip_obs(records: vec_deque::Iter<bam::Record>) {

    }
}

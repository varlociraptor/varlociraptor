//! Mismatches of reads against the reference, recorded per observation for the models of calling.
//!
//! A read from a source the reference lacks (a paralog, pseudogene or contaminant) shows mismatches
//! at about the divergence of its source, while a read from the mapped locus deviates only by its
//! base errors and by the variation that the other reads of the locus share.
//! Preprocessing records, per read, each mismatch with its base quality and how many of the other
//! fragments covering the position show the same base, so that calling can weigh them with its prior knowledge.

use std::collections::HashMap;
use std::fmt;
use std::ops::Range;
use std::str::FromStr;
use std::sync::Arc;

use anyhow::{anyhow, Result};
use rust_htslib::bam::{self, ext::BamRecordExtensions};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::utils;
use crate::variants::evidence::observations::read_observation::Evidence;

/// Mismatches stored per read; those of the highest quality come first.
pub const MAX_STORED_MISMATCHES: usize = 32;

/// The mismatches of a fragment against the reference.
/// Serialized as the string `bases;more;quality:shared/covering,...`, which suits both the BCF and the CSV output.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Mismatches {
    /// Aligned bases compared with the reference.
    pub bases: u32,
    /// Per mismatch, the base quality (PHRED), the number of other fragments showing the same base,
    /// and the number of other fragments covering the position;
    /// at most MAX_STORED_MISMATCHES entries, the highest qualities first.
    pub mismatches: Vec<(u8, u16, u16)>,
    /// Mismatches beyond the stored ones.
    pub more: u32,
}

impl Mismatches {
    /// All mismatches, stored or not.
    pub fn len(&self) -> usize {
        self.mismatches.len() + self.more as usize
    }
}

impl fmt::Display for Mismatches {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{};{};", self.bases, self.more)?;
        for (i, (qual, shared, covering)) in self.mismatches.iter().enumerate() {
            if i > 0 {
                write!(f, ",")?;
            }
            write!(f, "{qual}:{shared}/{covering}")?;
        }
        Ok(())
    }
}

impl FromStr for Mismatches {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let invalid = || anyhow!("invalid mismatch summary: {s}");
        let mut fields = s.split(';');
        let bases = fields
            .next()
            .and_then(|v| v.parse().ok())
            .ok_or_else(invalid)?;
        let more = fields
            .next()
            .and_then(|v| v.parse().ok())
            .ok_or_else(invalid)?;
        let mismatches = fields
            .next()
            .ok_or_else(invalid)?
            .split(',')
            .filter(|item| !item.is_empty())
            .map(|item| {
                let (qual, counts) = item.split_once(':').ok_or_else(invalid)?;
                let (shared, covering) = counts.split_once('/').ok_or_else(invalid)?;
                Ok((
                    qual.parse().map_err(|_| invalid())?,
                    shared.parse().map_err(|_| invalid())?,
                    covering.parse().map_err(|_| invalid())?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Mismatches {
            bases,
            mismatches,
            more,
        })
    }
}

impl Serialize for Mismatches {
    fn serialize<S: Serializer>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for Mismatches {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> std::result::Result<Self, D::Error> {
        String::deserialize(deserializer)?
            .parse()
            .map_err(serde::de::Error::custom)
    }
}

/// Whether the record is a consensus read: duplex (`ac`, `bc`) or single strand (`cD`) consensus of fgbio or fgumi,
/// or a consensus read of rust-bio-tools (`SI`).
pub(crate) fn is_consensus_read(record: &bam::Record) -> bool {
    record.aux(b"ac").is_ok()
        || record.aux(b"bc").is_ok()
        || record.aux(b"cD").is_ok()
        || utils::aux_tag_strand_info(record).is_some()
}

fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// An observed base of a fragment: its index (A, C, G, T), quality, and whether it deviates from the reference.
#[derive(Clone, Copy)]
struct FragmentBase {
    base: usize,
    qual: u8,
    mismatch: bool,
}

/// The bases of a fragment by reference position, each marked as a match or a mismatch against the reference.
/// Bases that are N in the read or in the reference are skipped.
/// Both mates of a pair read the same molecule, so a position covered by both counts once, with the better quality.
fn fragment_bases(evidence: &Evidence, ref_seq: &[u8]) -> HashMap<u64, FragmentBase> {
    let mut bases = HashMap::new();
    for record in evidence.records() {
        let seq = record.seq();
        let qual = record.qual();
        for [qpos, rpos] in record.aligned_pairs() {
            let (qpos, rpos) = (qpos as usize, rpos as usize);
            let ref_base = match ref_seq.get(rpos).map(|base| base.to_ascii_uppercase()) {
                Some(ref_base) if base_index(ref_base).is_some() => ref_base,
                _ => continue,
            };
            let read_base = seq[qpos].to_ascii_uppercase();
            if let Some(base) = base_index(read_base) {
                let observed = FragmentBase {
                    base,
                    qual: qual[qpos],
                    mismatch: read_base != ref_base,
                };
                let entry = bases.entry(rpos as u64).or_insert(observed);
                if observed.qual > entry.qual {
                    *entry = observed;
                }
            }
        }
    }
    bases
}

/// Per locus, how many fragments show each base at each position.
#[derive(Debug)]
pub(crate) struct MismatchCounts {
    /// Sequence of the contig of the locus.
    ref_seq: Arc<Vec<u8>>,
    /// Per reference position, the number of fragments covering it (index 4) and showing each base (indices 0 to 3).
    counts: HashMap<u64, [u32; 5]>,
    /// Reference ranges whose bases are left to the allele likelihoods.
    excluded: Vec<Range<u64>>,
}

impl MismatchCounts {
    pub(crate) fn new(
        evidences: &[Evidence],
        ref_seq: Arc<Vec<u8>>,
        excluded: Vec<Range<u64>>,
    ) -> Self {
        let mut counts: HashMap<u64, [u32; 5]> = HashMap::new();
        for evidence in evidences {
            for (rpos, observed) in fragment_bases(evidence, &ref_seq) {
                let entry = counts.entry(rpos).or_insert([0; 5]);
                entry[observed.base] += 1;
                entry[4] += 1;
            }
        }
        MismatchCounts {
            ref_seq,
            counts,
            excluded,
        }
    }

    fn excluded(&self, rpos: u64) -> bool {
        self.excluded.iter().any(|range| range.contains(&rpos))
    }

    /// The mismatches of the evidence outside the excluded ranges.
    pub(crate) fn summary(&self, evidence: &Evidence) -> Mismatches {
        let mut bases = 0u32;
        let mut mismatches = Vec::new();
        for (rpos, observed) in fragment_bases(evidence, &self.ref_seq) {
            if self.excluded(rpos) {
                continue;
            }
            bases += 1;
            if observed.mismatch {
                let counts = &self.counts[&rpos];
                let count = |n: u32| n.saturating_sub(1).min(u16::MAX as u32) as u16;
                mismatches.push((
                    observed.qual,
                    count(counts[observed.base]),
                    count(counts[4]),
                ));
            }
        }
        mismatches.sort_unstable_by(|a, b| b.cmp(a));
        let more = mismatches.len().saturating_sub(MAX_STORED_MISMATCHES) as u32;
        mismatches.truncate(MAX_STORED_MISMATCHES);
        Mismatches {
            bases,
            mismatches,
            more,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::evidence::observations::read_observation::AlignmentRecord;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn record(cigar: Vec<Cigar>, seq: &[u8], qual: u8) -> bam::Record {
        let mut record = bam::Record::new();
        let cigar = CigarString(cigar);
        record.set(b"read", Some(&cigar), seq, &vec![qual; seq.len()]);
        record.set_pos(2);
        record
    }

    fn evidence(record: bam::Record) -> Evidence {
        Evidence::SingleEndSequencingRead(AlignmentRecord::new(Arc::new(record), None))
    }

    #[test]
    fn test_fragment_bases() {
        // reference from position 2: ACGTAAGGTCCAN; query: 3 soft clipped, 4 matched (mismatch at the third),
        // deletion of 2, 5 matched (mismatch at the first), insertion, 3 matched
        // (the second against the N of the reference, the third past its end).
        let ref_seq = b"NNACGTAAGGTCCAN";
        let rec = record(
            vec![
                Cigar::SoftClip(3),
                Cigar::Match(4),
                Cigar::Del(2),
                Cigar::Match(5),
                Cigar::Ins(1),
                Cigar::Match(3),
            ],
            b"GGGACATTGTCCAACG",
            40,
        );
        let bases = fragment_bases(&evidence(rec), ref_seq);
        let mismatches: Vec<u64> = (2..15)
            .filter(|rpos| bases.get(rpos).is_some_and(|base| base.mismatch))
            .collect();
        assert_eq!(mismatches, vec![4, 8]);
        // the deleted positions, the N and the position past the reference end are not observed
        assert!(!bases.contains_key(&6) && !bases.contains_key(&7) && !bases.contains_key(&14));
        assert_eq!(bases.len(), 10);
    }

    #[test]
    fn test_summary() {
        let ref_seq: Arc<Vec<u8>> = Arc::new(b"NNACGTACGTACGT".to_vec());
        // three reads over positions 2 to 13: two share a mismatch at 5, one has another at 9
        let reads = vec![
            evidence(record(vec![Cigar::Match(12)], b"ACGAACGTACGT", 30)),
            evidence(record(vec![Cigar::Match(12)], b"ACGAACGTACGT", 30)),
            evidence(record(vec![Cigar::Match(12)], b"ACGTACGCACGT", 20)),
        ];
        let counts = MismatchCounts::new(&reads, Arc::clone(&ref_seq), vec![]);
        let first = counts.summary(&reads[0]);
        assert_eq!(first.bases, 12);
        assert_eq!(first.mismatches, vec![(30, 1, 2)]);
        assert_eq!(first.more, 0);
        let third = counts.summary(&reads[2]);
        assert_eq!(third.mismatches, vec![(20, 0, 2)]);
        // the excluded range hides the mismatch and its bases
        let counts = MismatchCounts::new(&reads, ref_seq, vec![4..7]);
        let first = counts.summary(&reads[0]);
        assert_eq!(first.bases, 9);
        assert!(first.mismatches.is_empty());
    }

    #[test]
    fn test_mismatches_string() {
        let mismatches = Mismatches {
            bases: 150,
            mismatches: vec![(90, 0, 180), (12, 7, 175)],
            more: 3,
        };
        let text = mismatches.to_string();
        assert_eq!(text, "150;3;90:0/180,12:7/175");
        assert_eq!(text.parse::<Mismatches>().unwrap(), mismatches);
        assert_eq!(
            "12;0;".parse::<Mismatches>().unwrap(),
            Mismatches {
                bases: 12,
                ..Default::default()
            }
        );
        assert!("12;0".parse::<Mismatches>().is_err());
    }
}

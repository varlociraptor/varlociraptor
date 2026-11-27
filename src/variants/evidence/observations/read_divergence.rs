//! Mapping of reads from their own mismatches.
//!
//! MAPQ only reflects alternative placements within the reference.
//! A read may also come from a source the reference lacks
//! (a paralog, pseudogene or contaminant, sequence that is absent or divergent in the reference).
//! Such a read shows mismatches against the reference at about the divergence of its source,
//! while a read from the mapped locus deviates only by its base errors
//! and by the variation that the other reads of the locus share.
//! For consensus reads, whose base errors are far below any divergence,
//! this makes the read's own mismatches the stronger mapping evidence;
//! for raw reads the base qualities keep the evidence weak,
//! which is why the model is applied to consensus reads by default.

use std::collections::HashMap;
use std::ops::Range;
use std::sync::Arc;

use bio::stats::LogProb;
use rust_htslib::bam::{self, ext::BamRecordExtensions};

use crate::variants::evidence::observations::read_observation::{AlignmentRecord, Evidence};

/// Reads to which the divergence model applies.
#[derive(
    Display,
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    EnumString,
    EnumIter,
    IntoStaticStr,
    VariantNames,
    Default,
)]
#[strum(serialize_all = "kebab_case")]
#[serde(rename_all = "kebab-case")]
pub enum ReadDivergenceReads {
    /// Consensus reads, recognized by the tags `ac`, `bc` or `cD` as written by fgbio or fgumi.
    #[default]
    Auto,

    /// All reads.
    Always,

    /// No reads.
    Never,
}

/// Settings of the divergence model, per sample.
#[derive(Debug, Clone, Copy)]
pub struct ReadDivergenceSettings {
    pub reads: ReadDivergenceReads,

    /// Divergence of the source a mismapped read would come from.
    pub divergence: f64,

    /// A base shared by at least this fraction of the other fragments of the locus is variation of the sample;
    /// a rarer shared base is divergence,
    /// since a rare haplotype with several linked differences is far less likely than reads from another source,
    /// while germline variation sits at one over the ploidy and above.
    pub variation_min_frequency: f64,

    /// Bases this far around a candidate that is realigned (anything but SNVs and MNVs) are left out,
    /// since their alignment follows the aligner's placement of the candidate.
    pub exclusion_window: u64,
}

impl Default for ReadDivergenceSettings {
    fn default() -> Self {
        ReadDivergenceSettings {
            reads: ReadDivergenceReads::Auto,
            divergence: 0.03,
            variation_min_frequency: 0.2,
            exclusion_window: 64,
        }
    }
}

/// Whether the record carries consensus tags: `ac` and `bc` of a duplex read or `cD` of a single-strand consensus read,
/// as written by fgbio or fgumi.
pub(crate) fn is_consensus_read(record: &bam::Record) -> bool {
    record.aux(b"ac").is_ok() || record.aux(b"bc").is_ok() || record.aux(b"cD").is_ok()
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

/// An observed base of a fragment: its index (A, C, G, T), error probability,
/// and whether it deviates from the reference.
#[derive(Clone, Copy)]
struct FragmentBase {
    base: usize,
    prob_error: f64,
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
                    prob_error: 10f64.powf(-(qual[qpos] as f64) / 10.0),
                    mismatch: read_base != ref_base,
                };
                let entry = bases.entry(rpos as u64).or_insert(observed);
                if observed.prob_error < entry.prob_error {
                    *entry = observed;
                }
            }
        }
    }
    bases
}

/// Per-locus model of read divergence: which mismatches the mapped locus explains, and the rate of the rest.
#[derive(Debug)]
pub(crate) struct ReadDivergenceModel {
    settings: ReadDivergenceSettings,

    /// Sequence of the contig of the locus.
    ref_seq: Arc<Vec<u8>>,

    /// Per reference position, the number of fragments covering it (index 4) and showing each base (indices 0 to 3).
    counts: HashMap<u64, [u32; 5]>,

    /// Rate of mismatches per aligned base that no other fragment of the locus shares.
    prob_background: f64,

    /// Reference ranges whose bases are left to the allele likelihoods.
    excluded: Vec<Range<u64>>,
}

impl ReadDivergenceModel {
    fn applies(&self, record: &bam::Record) -> bool {
        match self.settings.reads {
            ReadDivergenceReads::Auto => is_consensus_read(record),
            ReadDivergenceReads::Always => true,
            ReadDivergenceReads::Never => false,
        }
    }

    fn excluded(&self, rpos: u64) -> bool {
        self.excluded.iter().any(|range| range.contains(&rpos))
    }

    /// Build the model from all evidence of a locus, or None if the model applies to none of it.
    pub(crate) fn new(
        evidences: &[Evidence],
        ref_seq: Arc<Vec<u8>>,
        excluded: Vec<Range<u64>>,
        settings: ReadDivergenceSettings,
    ) -> Option<Self> {
        let mut model = ReadDivergenceModel {
            settings,
            ref_seq,
            counts: HashMap::new(),
            prob_background: 0.0,
            excluded,
        };
        let mut fragments = Vec::new();
        let mut any = false;
        for evidence in evidences {
            if !evidence.records().all(|record| model.applies(record)) {
                continue;
            }
            any = true;
            let bases = fragment_bases(evidence, &model.ref_seq);
            for (rpos, observed) in &bases {
                let counts = model.counts.entry(*rpos).or_insert([0; 5]);
                counts[observed.base] += 1;
                counts[4] += 1;
            }
            fragments.push(bases);
        }
        if !any {
            return None;
        }
        // METHOD: the background is the rate of mismatches that are not variation of the sample,
        // with one pseudo observation each (a flat prior on the rate).
        let mut n_bases = 1usize;
        let mut n_unexplained = 1usize;
        for bases in &fragments {
            for (rpos, observed) in bases {
                if model.excluded(*rpos) {
                    continue;
                }
                n_bases += 1;
                if observed.mismatch {
                    let counts = &model.counts[rpos];
                    let shared = if counts[4] > 1 {
                        (counts[observed.base] - 1) as f64 / (counts[4] - 1) as f64
                    } else {
                        0.0
                    };
                    if shared < settings.variation_min_frequency {
                        n_unexplained += 1;
                    }
                }
            }
        }
        model.prob_background = n_unexplained as f64 / n_bases as f64;
        Some(model)
    }

    /// Mapping probability of the evidence after updating the given one (from MAPQ)
    /// with the likelihood ratio of its bases under another source against the mapped locus,
    /// and the probability of the observed allele at another source.
    ///
    /// Returns None if the model does not apply to the evidence.
    pub(crate) fn recalibrate(
        &self,
        evidence: &Evidence,
        prob_mapping: LogProb,
        alt_support: bool,
    ) -> Option<(LogProb, LogProb)> {
        if !evidence.records().all(|record| self.applies(record)) {
            return None;
        }
        let divergence = self.settings.divergence;
        let mut log_ratio = 0.0;
        for (rpos, observed) in fragment_bases(evidence, &self.ref_seq) {
            if self.excluded(rpos) {
                continue;
            }
            if observed.mismatch {
                // METHOD: the mapped locus explains the mismatch by a base error,
                // by the background of unshared mismatches,
                // or by the fraction of the other fragments that show the same base there
                // if that is variation of the sample (see variation_min_frequency);
                // another source by its divergence.
                let counts = &self.counts[&rpos];
                let shared = if counts[4] > 1 {
                    (counts[observed.base] - 1) as f64 / (counts[4] - 1) as f64
                } else {
                    0.0
                };
                let variation = if shared >= self.settings.variation_min_frequency {
                    shared
                } else {
                    0.0
                };
                let prob_mapped = observed.prob_error.max(self.prob_background).max(variation);
                log_ratio += (divergence / 3.0).ln() - prob_mapped.ln();
            } else {
                log_ratio += (1.0 - divergence).ln() - (1.0 - self.prob_background).ln();
            }
        }
        let prob_mismapping = prob_mapping.ln_one_minus_exp();
        let mismapped = LogProb(*prob_mismapping + log_ratio);
        let prob_mismapping = LogProb(*mismapped - *mismapped.ln_add_exp(prob_mapping));
        let prob_missed_allele = if alt_support {
            // METHOD: another source shows the alt base where it differs from the mapped locus.
            LogProb((divergence / 3.0).ln())
        } else {
            LogProb((1.0 - divergence).ln())
        };
        Some((prob_mismapping.ln_one_minus_exp(), prob_missed_allele))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn record(cigar: Vec<Cigar>, seq: &[u8]) -> bam::Record {
        let mut record = bam::Record::new();
        let cigar = CigarString(cigar);
        record.set(b"read", Some(&cigar), seq, &vec![40u8; seq.len()]);
        record.set_pos(2);
        record
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
        );
        let evidence = Evidence::SingleEndSequencingRead(AlignmentRecord::new(Arc::new(rec), None));
        let bases = fragment_bases(&evidence, ref_seq);
        let mismatches: Vec<u64> = (2..15)
            .filter(|rpos| bases.get(rpos).is_some_and(|base| base.mismatch))
            .collect();
        assert_eq!(mismatches, vec![4, 8]);
        // the deleted positions, the N and the position past the reference end are not observed
        assert!(!bases.contains_key(&6) && !bases.contains_key(&7) && !bases.contains_key(&14));
        assert_eq!(bases.len(), 10);
    }
}

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
use std::sync::atomic::{AtomicBool, Ordering};

use bio::stats::LogProb;
use rust_htslib::bam::{self, record::Aux, record::Cigar};

use crate::variants::evidence::observations::read_observation::Evidence;

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

static WARNED_MISSING_MD: AtomicBool = AtomicBool::new(false);

/// Whether the record carries consensus tags: `ac` and `bc` of a duplex read or `cD` of a single-strand consensus read,
/// as written by fgbio or fgumi.
pub(crate) fn is_consensus_read(record: &bam::Record) -> bool {
    record.aux(b"ac").is_ok() || record.aux(b"bc").is_ok() || record.aux(b"cD").is_ok()
}

/// Aligned bases of a record as (reference position, query position), via the CIGAR.
fn aligned_bases(record: &bam::Record) -> Vec<(u64, usize)> {
    let mut bases = Vec::with_capacity(record.seq_len());
    let mut rpos = record.pos() as u64;
    let mut qpos = 0;
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                for _ in 0..*n {
                    bases.push((rpos, qpos));
                    rpos += 1;
                    qpos += 1;
                }
            }
            Cigar::Ins(n) | Cigar::SoftClip(n) => qpos += *n as usize,
            Cigar::Del(n) | Cigar::RefSkip(n) => rpos += *n as u64,
            Cigar::HardClip(_) | Cigar::Pad(_) => (),
        }
    }
    bases
}

/// Query positions of the mismatches of a record against the reference, from its MD tag.
/// None if the record has no (or a malformed) MD tag.
fn mismatches(record: &bam::Record) -> Option<Vec<usize>> {
    let md = match record.aux(b"MD") {
        Ok(Aux::String(md)) => md.as_bytes(),
        _ => return None,
    };
    let mut i = 0;
    let number = |i: &mut usize| {
        let start = *i;
        while *i < md.len() && md[*i].is_ascii_digit() {
            *i += 1;
        }
        std::str::from_utf8(&md[start..*i])
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0)
    };
    let mut matches_left = number(&mut i);
    let mut mismatches = Vec::new();
    let mut qpos = 0;
    for op in record.cigar().iter() {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                for _ in 0..*n {
                    if matches_left == 0 {
                        if i < md.len() && md[i].is_ascii_alphabetic() {
                            mismatches.push(qpos);
                            i += 1;
                            matches_left = number(&mut i);
                        } else {
                            return None;
                        }
                    } else {
                        matches_left -= 1;
                    }
                    qpos += 1;
                }
            }
            Cigar::Ins(n) | Cigar::SoftClip(n) => qpos += *n as usize,
            Cigar::Del(_) => {
                if i < md.len() && md[i] == b'^' {
                    i += 1;
                    while i < md.len() && md[i].is_ascii_alphabetic() {
                        i += 1;
                    }
                    matches_left = number(&mut i);
                } else {
                    return None;
                }
            }
            Cigar::RefSkip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => (),
        }
    }
    Some(mismatches)
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

/// The bases of a fragment by reference position.
/// Both mates of a pair observe the same molecule,
/// so a position covered by both counts once, with the better quality.
fn fragment_bases(evidence: &Evidence) -> Option<HashMap<u64, FragmentBase>> {
    let mut bases = HashMap::new();
    for record in evidence.records() {
        let mismatches = mismatches(record)?;
        let seq = record.seq();
        let qual = record.qual();
        for (rpos, qpos) in aligned_bases(record) {
            if let Some(base) = base_index(seq[qpos]) {
                let observed = FragmentBase {
                    base,
                    prob_error: 10f64.powf(-(qual[qpos] as f64) / 10.0),
                    mismatch: mismatches.contains(&qpos),
                };
                let entry = bases.entry(rpos).or_insert(observed);
                if observed.prob_error < entry.prob_error {
                    *entry = observed;
                }
            }
        }
    }
    Some(bases)
}

/// Per-locus model of read divergence: which mismatches the mapped locus explains, and the rate of the rest.
#[derive(Debug)]
pub(crate) struct ReadDivergenceModel {
    settings: ReadDivergenceSettings,

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
        excluded: Vec<Range<u64>>,
        settings: ReadDivergenceSettings,
    ) -> Option<Self> {
        let mut model = ReadDivergenceModel {
            settings,
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
            match fragment_bases(evidence) {
                Some(bases) => {
                    for (rpos, observed) in &bases {
                        let counts = model.counts.entry(*rpos).or_insert([0; 5]);
                        counts[observed.base] += 1;
                        counts[4] += 1;
                    }
                    fragments.push(bases);
                }
                None => {
                    if !WARNED_MISSING_MD.swap(true, Ordering::Relaxed) {
                        warn!(
                            "Reads without MD tag found; the read divergence model needs MD tags \
                             (bwa writes them, minimap2 with --MD) and leaves such reads unchanged."
                        );
                    }
                }
            }
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
    /// Returns None if the model does not apply to the evidence or its records lack MD tags.
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
        for (rpos, observed) in fragment_bases(evidence)? {
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
    use rust_htslib::bam::record::{Aux, CigarString};

    fn record(cigar: Vec<Cigar>, seq: &[u8], md: &str) -> bam::Record {
        let mut record = bam::Record::new();
        let cigar = CigarString(cigar);
        record.set(b"read", Some(&cigar), seq, &vec![40u8; seq.len()]);
        record.set_pos(100);
        record.push_aux(b"MD", Aux::String(md)).unwrap();
        record
    }

    #[test]
    fn test_mismatches_from_md() {
        // query: 3 soft clipped, 4 matched (mismatch at the third), deletion,
        // 5 matched (mismatch at the first), insertion, 3 matched (mismatch at the first).
        let rec = record(
            vec![
                Cigar::SoftClip(3),
                Cigar::Match(4),
                Cigar::Del(2),
                Cigar::Match(5),
                Cigar::Ins(1),
                Cigar::Match(3),
            ],
            b"GGGACGTACGTAAGGG",
            "2A1^TT0G4A2",
        );
        assert_eq!(mismatches(&rec).unwrap(), vec![5, 7, 13]);
        assert!(mismatches(&record(vec![Cigar::Match(4)], b"ACGT", "2X")).is_none());
    }

    #[test]
    fn test_aligned_bases() {
        let bases = aligned_bases(&record(
            vec![
                Cigar::SoftClip(2),
                Cigar::Match(2),
                Cigar::Del(1),
                Cigar::Match(1),
            ],
            b"GGACT",
            "3",
        ));
        assert_eq!(bases, vec![(100, 2), (101, 3), (103, 4)]);
    }
}

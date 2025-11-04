use super::ToVariantRepresentation;
use crate::variants::evidence::realignment::Realignable;
use crate::variants::model;
use crate::{estimation::alignment_properties::AlignmentProperties, variants::sample::Readtype};

use super::MultiLocus;
use crate::variants::evidence::bases::prob_read_base;
use crate::variants::evidence::observations::read_observation::{AlignmentRecord, Strand};
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, Evidence, Overlap, SingleLocus, Variant,
};
use anyhow::Result;
use bio::stats::{LogProb, Prob};
use bio_types::genome::{self, AbstractInterval, AbstractLocus};
use log::warn;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;
use std::collections::HashMap;
use std::rc::Rc;

#[derive(Debug)]
pub(crate) struct Methylation {
    loci: MultiLocus,
    readtype: Readtype,
}

impl Methylation {
    pub(crate) fn new(locus: genome::Locus, readtype: Readtype) -> Self {
        Methylation {
            loci: MultiLocus::from_single_locus(SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + 1,
            ))),
            readtype,
        }
    }

    fn allele_support_per_read(
        &self,
        read: &AlignmentRecord,
        is_long_read: bool,
    ) -> Result<Option<AlleleSupport>> {
        let mut position = self.locus().range().start;
        if read_reverse_orientation(read) {
            position += 1;
        }
        if let Some(qpos) = read
            .cigar_cached()
            .unwrap()
            // TODO expect u64 in read_pos
            .read_pos(position as u32, false, false)?
        {
            if let Some((prob_alt, prob_ref)) =
                process_read(read, read.prob_methylation(), qpos, is_long_read)
            {
                let strand = if prob_ref != prob_alt {
                    Strand::from_record_and_pos(read, qpos as usize)?
                } else {
                    // METHOD: if record is not informative, we don't want to
                    // retain its information (e.g. strand).
                    Strand::no_strand_info()
                };
                Ok(Some(
                    AlleleSupportBuilder::default()
                        .prob_ref_allele(prob_ref)
                        .prob_alt_allele(prob_alt)
                        .strand(strand)
                        .read_position(Some(qpos))
                        .third_allele_evidence(None)
                        .build()
                        .unwrap(),
                ))
            } else {
                Ok(None)
            }
        } else {
            // a read that shows methylation might have the respective position in the
            // reference skipped (Cigar op 'N'), and the library should not choke on those reads
            // but instead needs to know NOT to add those reads (as observations) further up
            Ok(None)
        }
    }

    fn locus(&self) -> &SingleLocus {
        &self.loci[0]
    }
}

/// Checks if the read has MM information
///
/// # Returns
///
/// bool: Read has MM information
fn mm_exists(evidence: &Evidence) -> bool {
    match evidence {
        Evidence::SingleEndSequencingRead(record) => {
            // Check for MM tags in the single-end sequencing read
            mm_tag_exists(record)
        }
        Evidence::PairedEndSequencingRead { left, right } => {
            // Check for MM tags in both left and right reads
            mm_tag_exists(left) || mm_tag_exists(right)
        }
    }
}

/// Helper function to check MM tags for a single record
fn mm_tag_exists(record: &Rc<bam::Record>) -> bool {
    matches!(
        (record.aux(b"Mm"), record.aux(b"MM")),
        (Ok(_), _) | (_, Ok(_))
    )
}

fn is_5mc_header(header: &str) -> bool {
    header.starts_with("C+m") || header.starts_with("C-m")
}
/// Computes the positions and probabilities of methylated bases in PacBio and Nanopore read data.
/// Handles multiple MM blocks (e.g. A+a., C+h., etc.)
///
/// # Returns
/// pos_methylated_bases: Vector of positions (0-based read indices) of methylated bases
///
///
pub fn extract_mm_ml_5mc(read: &Rc<Record>) -> Option<HashMap<usize, LogProb>> {
    let mm_tag = match (read.aux(b"Mm"), read.aux(b"MM")) {
        (Ok(tag), _) => tag,
        (_, Ok(tag)) => tag,
        _ => return None,
    };

    let ml_tag = match (read.aux(b"Ml"), read.aux(b"ML")) {
        (Ok(tag), _) => tag,
        (_, Ok(tag)) => tag,
        _ => return None,
    };

    let mm = if let Aux::String(tag_value) = mm_tag {
        tag_value.to_owned()
    } else {
        return None;
    };

    let ml = if let Aux::ArrayU8(tag_value) = ml_tag {
        tag_value
    } else {
        return None;
    };

    let read_seq = read.seq().as_bytes();

    let mut pos_to_prob: HashMap<usize, LogProb> = HashMap::new();
    let mut ml_index = 0;

    for block in mm.split(';') {
        if block.is_empty() {
            continue;
        }

        let (header, positions_str) = match block.split_once(',') {
            Some(x) => x,
            None => continue,
        };

        // Modified bases in the read
        // The MM tag encodes positions as deltas from the previous position
        let methylated_bases: Vec<usize> = positions_str
            .split(',')
            .filter_map(|s| s.parse::<usize>().ok())
            .collect();

        if is_5mc_header(header) {
            // Positions of 'C' bases in the read
            let mut pos_read_base: Vec<usize> = read_seq
                .iter()
                .enumerate()
                .filter(|&(_, &c)| {
                    (c == b'C' && !read_reverse_orientation(read))
                        || (c == complement_base(b'C') && read_reverse_orientation(read))
                })
                .map(|(i, _)| i)
                .collect();

            if read_reverse_orientation(read) {
                pos_read_base.reverse();
            }

            let mut meth_pos = 0;
            for next_base in methylated_bases {
                // Specific position of methylated base in the read
                meth_pos += next_base;
                // In a few cases, the MM tag might contain positions that exceed the read length
                if meth_pos <= pos_read_base.len() {
                    let abs_pos = pos_read_base[meth_pos];
                    let prob_val = ml.get(ml_index).unwrap_or(0);
                    let prob = LogProb::from(Prob((f64::from(prob_val) + 0.5) / 256.0));

                    pos_to_prob.insert(abs_pos, prob);
                }
                ml_index += 1;
                meth_pos += 1; // Move to the next position
            }
        } else {
            // Skip other modification types (no 5mC)
            ml_index += methylated_bases.len();
        }
    }
    let mut items: Vec<_> = pos_to_prob.iter().collect();

    // Sort by position in descending order
    items.sort_by(|a, b| b.0.cmp(a.0));
    dbg!(&items);
    Some(pos_to_prob)
}

// Returns the complement base for a given base
fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => base,
    }
}

fn process_read(
    read: &Rc<Record>,
    meth_info: &Option<Rc<HashMap<usize, LogProb>>>,
    qpos: u32,
    is_long_read: bool,
) -> Option<(LogProb, LogProb)> {
    if (is_long_read && meth_info.is_none())
        || mutation_occurred(read_reverse_orientation(read), read, qpos, is_long_read)
        || read_invalid(read.inner.core.flag)
    {
        return None;
    }

    if is_long_read {
        Some(compute_probs_long_read(meth_info.as_ref().unwrap(), qpos))
    } else {
        Some(compute_probs_short_read(
            read_reverse_orientation(read),
            read,
            qpos,
        ))
    }
}

/// Computes the probability of methylation/no methylation of a given position in an PacBio/ Nanopore read. Takes mapping probability into account
///
/// # Returns
///
/// prob_alt: Probability of methylation (alternative)
/// prob_ref: Probability of no methylation (reference)
fn compute_probs_long_read(
    pos_to_probs: &HashMap<usize, LogProb>,
    qpos: u32,
) -> (LogProb, LogProb) {
    // let pos_in_read = qpos ;
    let prob_alt;
    let prob_ref;
    if let Some(value) = pos_to_probs.get(&(qpos as usize)) {
        prob_alt = value.to_owned();
        prob_ref = LogProb::from(Prob(1_f64 - prob_alt.0.exp()));
    } else {
        prob_alt = LogProb::from(Prob(0.0));
        prob_ref = LogProb::from(Prob(1.0));
    }
    (prob_alt, prob_ref)
}

/// Computes the probability of methylation/no methylation of a given position in an Illumina read. Takes mapping probability into account
///
/// # Returns
///
/// prob_alt: Probability of methylation (alternative)
/// prob_ref: Probability of no methylation (reference)
pub fn compute_probs_short_read(
    read_reverse: bool,
    record: &Rc<Record>,
    qpos: u32,
) -> (LogProb, LogProb) {
    let (ref_base, bisulfite_base) = if !read_reverse {
        (b'C', b'T')
    } else {
        (b'G', b'A')
    };
    let read_base = unsafe { record.seq().decoded_base_unchecked(qpos as usize) };
    let base_qual = unsafe { *record.qual().get_unchecked(qpos as usize) };
    let prob_alt = prob_read_base(read_base, ref_base, base_qual);
    let prob_ref = prob_read_base(read_base, bisulfite_base, base_qual);
    (prob_alt, prob_ref)
}

fn mutation_occurred(
    read_reverse: bool,
    record: &Rc<Record>,
    qpos: u32,
    is_long_read: bool,
) -> bool {
    let (read_base, mutation_bases) = if read_reverse {
        let read_base = unsafe { record.seq().decoded_base_unchecked(qpos as usize) };
        let mutation_bases = if is_long_read {
            vec![b'C', b'A', b'T']
        } else {
            vec![b'C', b'T']
        };
        (read_base, mutation_bases)
    } else {
        let read_base = unsafe { record.seq().decoded_base_unchecked(qpos as usize) };
        let mutation_bases = if is_long_read {
            vec![b'G', b'A', b'T']
        } else {
            vec![b'A', b'G']
        };
        (read_base, mutation_bases)
    };

    if mutation_bases.contains(&read_base) {
        warn!(
            "The record {:?} on position {:?} is not considered because a mutation occurred",
            String::from_utf8_lossy(record.qname()),
            qpos
        );
        return true;
    }
    false
}

/// Computes if a given read is valid (Right now we only accept specific flags)
///
/// # Returns
///
/// bool: True, if read is valid, else false
fn read_invalid(flag: u16) -> bool {
    !matches!(flag, 0 | 16 | 83 | 99 | 147 | 163)
}

impl Variant for Methylation {
    fn is_imprecise(&self) -> bool {
        false
    }

    fn loci(&self) -> &MultiLocus {
        &self.loci
    }

    /// Determine whether the evidence is suitable to assessing probabilities
    /// (i.e. overlaps the variant in the right way).
    ///
    /// # Returns
    ///
    /// The index of the loci for which this evidence is valid, `None` if invalid.
    fn is_valid_evidence(
        &self,
        evidence: &Evidence,
        _: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        if match evidence {
            Evidence::SingleEndSequencingRead(read) => {
                let start_offset = if read_reverse_orientation(read) { 1 } else { 0 };
                matches!(
                    self.locus().overlap(read, false, start_offset, 0),
                    Overlap::Enclosing
                )
            }
            Evidence::PairedEndSequencingRead { left, right } => {
                let start_offset_left = if read_reverse_orientation(left) { 1 } else { 0 };
                let start_offset_right = if read_reverse_orientation(right) {
                    1
                } else {
                    0
                };

                matches!(
                    self.locus().overlap(left, false, start_offset_left, 0),
                    Overlap::Enclosing
                ) || matches!(
                    self.locus().overlap(right, false, start_offset_right, 0),
                    Overlap::Enclosing
                )
            }
        } {
            if match self.readtype {
                // Some single PacBio reads don't have an MM:Z value and are therefore not legal
                Readtype::Illumina => true,
                Readtype::PacBio => mm_exists(evidence),
                Readtype::Nanopore => mm_exists(evidence),
            } {
                Some(vec![0])
            } else {
                None
            }
        } else {
            None
        }
    }

    fn allele_support(
        &self,
        evidence: &Evidence,
        _alignment_properties: &AlignmentProperties,
        _alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        match evidence {
            Evidence::SingleEndSequencingRead(read) => match self.readtype {
                Readtype::Illumina => Ok(self.allele_support_per_read(read, false)?),
                Readtype::PacBio | Readtype::Nanopore => {
                    Ok(self.allele_support_per_read(read, true)?)
                }
            },
            Evidence::PairedEndSequencingRead { left, right } => match self.readtype {
                Readtype::Illumina => {
                    let left_support = self.allele_support_per_read(left, false)?;
                    let right_support = self.allele_support_per_read(right, false)?;

                    match (left_support, right_support) {
                        (Some(mut left_support), Some(right_support)) => {
                            left_support.merge(&right_support);
                            Ok(Some(left_support))
                        }
                        (Some(left_support), None) => Ok(Some(left_support)),
                        (None, Some(right_support)) => Ok(Some(right_support)),
                        (None, None) => Ok(None),
                    }
                }
                Readtype::PacBio | Readtype::Nanopore => {
                    let left_support = self.allele_support_per_read(left, true)?;
                    let right_support = self.allele_support_per_read(right, true)?;

                    match (left_support, right_support) {
                        (Some(mut left_support), Some(right_support)) => {
                            left_support.merge(&right_support);
                            Ok(Some(left_support))
                        }
                        (Some(left_support), None) => Ok(Some(left_support)),
                        (None, Some(right_support)) => Ok(Some(right_support)),
                        (None, None) => Ok(None),
                    }
                }
            },
        }
    }

    /// Calculate probability to sample a record length like the given one from the alt allele.
    fn prob_sample_alt(&self, _: &Evidence, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

impl ToVariantRepresentation for Methylation {
    fn to_variant_representation(&self) -> model::Variant {
        model::Variant::Methylation()
    }
}

pub(crate) fn read_reverse_orientation(read: &Rc<Record>) -> bool {
    let read_paired = read.is_paired();
    let read_reverse = read.is_reverse();
    let read_first = read.is_first_in_template();
    if read_paired {
        read_reverse && read_first || !read_reverse && !read_first
    } else {
        read_reverse
    }
}

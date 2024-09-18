use super::ToVariantRepresentation;
use crate::variants::evidence::realignment::Realignable;
use crate::variants::model;
use crate::{estimation::alignment_properties::AlignmentProperties, variants::sample::Readtype};

use crate::variants::evidence::bases::prob_read_base;
use crate::variants::evidence::observations::read_observation::Strand;
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, PairedEndEvidence, SingleEndEvidence, SingleLocus, Variant,
};
use rust_htslib::bam::Record;
use std::collections::HashMap;
use std::rc::Rc;

use anyhow::Result;
use bio::stats::{LogProb, Prob};
use bio_types::genome::{self, AbstractInterval, AbstractLocus};
use lazy_static::lazy_static;
use log::warn;
use rust_htslib::bam::record::Aux;
use std::sync::Mutex;

// We save methylation info of the single reads for PacBio and Nanopore in order to not recompute the information for every candidate
lazy_static! {
    static ref READ_TO_METH_PROBS: Mutex<HashMap<String, HashMap<usize, f64>>> =
        Mutex::new(HashMap::new());
}

#[derive(Debug)]
pub(crate) struct Methylation {
    locus: SingleLocus,
    readtype: Readtype,
}

impl Methylation {
    pub(crate) fn new(locus: genome::Locus, readtype: Readtype) -> Self {
        Methylation {
            locus: SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + 2,
            )),
            readtype,
        }
    }
}

/// Looks if the read has MM information
///
/// # Returns
///
/// bool: Read has MM information
fn mm_exist(read: &SingleEndEvidence) -> bool {
    let mm_tag_exists = match (read.aux(b"Mm"), read.aux(b"MM")) {
        (Ok(_), _) | (_, Ok(_)) => true, // True, wenn einer der Tags existiert
        _ => false,                      // False, wenn keiner der Tags existiert
    };
    mm_tag_exists
}

/// Computes the positions of methylated Cs in PacBio and Nanopore read data
///
/// # Returns
///
/// pos_methylated_cs: Vector of positions of methylated Cs in Read
pub fn meth_pos(read: &Rc<Record>) -> Option<Vec<usize>> {
    let mm_tag = match (read.aux(b"Mm"), read.aux(b"MM")) {
        (Ok(tag), _) => tag,
        (_, Ok(tag)) => tag,
        _ => {
            warn!(
                "MM value not found on chrom {:?}, pos {:?}",
                read.inner.core.tid, read
            );
            return None;
        }
    };

    if let Aux::String(tag_value) = mm_tag {
        let mut mm = tag_value.to_owned();

        if !mm.is_empty() {
            let read_reverse = SingleLocus::read_reverse_strand(read.inner.core.flag);
            let read_seq = String::from_utf8_lossy(&read.seq().as_bytes()).to_string();

            // Compute the positions of all Cs (or Gs for reverse strand) in the read
            // With MM tags it is important to know whether it is a forward or reverse strand. With forward, you go from front to back and count the Cs. With reverse, you go from back to front in the sequence, in the MM tag from front to back and count the Gs
            let mut pos_read_base: Vec<usize> = read_seq
                .chars()
                .enumerate()
                .filter(|&(_, c)| (c == 'C' && !read_reverse) || (c == 'G' && read_reverse))
                .map(|(index, _)| index)
                .collect();

            if read_reverse {
                pos_read_base.reverse()
            }

            // Compute which Cs are methylated with help of the MM-tag
            mm.pop();
            if let Some(remaining) = mm.strip_prefix("C+m") {
                // There is no convention if the MM tag starts with C+m, or C+m.,
                let methylated_part = if remaining.starts_with(',') {
                    remaining.strip_prefix(',').unwrap_or("")
                } else {
                    remaining.strip_prefix(".,").unwrap_or("")
                };

                let mut meth_pos = 0;

                let mut pos_methylated_cs: Vec<usize> = Vec::new(); // Initialize empty vector

                for position_str in methylated_part.split(',') {
                    let position: usize = match position_str.parse::<usize>() {
                        Ok(position) => position,
                        Err(_) => {
                            warn!("Invalid position format: {}", position_str);
                            return None; // Return empty vector on invalid format
                        }
                    };

                    meth_pos += position + 1;

                    if meth_pos <= pos_read_base.len() {
                        pos_methylated_cs.push(pos_read_base[meth_pos - 1]);
                    } else {
                        warn!(
                            "MM tag too long for read on id {:?}, chrom {:?}, pos {:?}",
                            String::from_utf8_lossy(read.qname()),
                            read.inner.core.tid,
                            read.inner.core.pos
                        );
                        return None; // Return empty vector on out-of-bounds index
                    }
                }

                if read_reverse {
                    pos_methylated_cs.reverse();
                }
                return Some(pos_methylated_cs);
            }
        }
        // No methylation info in this read
        else {
            return Some(Vec::new());
        }
    } else {
        warn!(
            "MM tag in bam file is not valid on chrom {:?}, pos {:?}",
            read.inner.core.tid, read.inner.core.pos
        );
        return None;
    }
    warn!(
        "Error while obtaining MM:Z tag on chrom {:?}, pos {:?}",
        read.inner.core.tid, read.inner.core.pos
    );
    None
}

/// Computes the probability of methylation in PacBio and Nanopore read data
///
/// # Returns
///
/// ml: Vector of methylation probabilities
pub fn meth_probs(read: &Rc<Record>) -> Option<Vec<LogProb>> {
    let ml_tag = match (read.aux(b"Ml"), read.aux(b"ML")) {
        (Ok(tag), _) => tag,
        (_, Ok(tag)) => tag,
        _ => {
            warn!(
                "ML value not found on chrom {:?}, pos {:?}",
                read.inner.core.tid, read.inner.core.pos
            );
            return None;
        }
    };
    let read_reverse = SingleLocus::read_reverse_strand(read.inner.core.flag);
    if let Aux::ArrayU8(tag_value) = ml_tag {
        let mut ml: Vec<LogProb> = tag_value
            .iter()
            .map(|val| LogProb::from(Prob((f64::from(val) + 0.5) / 256.0)))
            .collect();
        if read_reverse {
            ml.reverse();
        }
        Some(ml)
    } else {
        warn!(
            "MM tag in bam file is not valid on chrom {:?}, pos {:?}",
            read.inner.core.tid, read.inner.core.pos
        );
        None
    }
}

/// Position of CpG site in the read
///
/// # Returns
///
/// Option<position>: Position of CpG site if it is included in read
/// None else
fn get_qpos(read: &Rc<Record>, locus: &SingleLocus) -> Option<i32> {
    read.cigar_cached()
        .unwrap()
        .read_pos(locus.range().start as u32, false, false)
        .unwrap()
        .map(|qpos| qpos as i32)
}

fn process_read_illumina(read: &Rc<Record>, locus: &SingleLocus) -> Option<(LogProb, LogProb)> {
    let qpos = get_qpos(read, locus)?;
    let read_reverse = SingleLocus::read_reverse_strand(read.inner.core.flag);
    let c_not_included = qpos as usize == read.seq().len() - 1 && read_reverse;

    let is_invalid = read_invalid(read.inner.core.flag);
    let mutation_occurred = mutation_occurred_illumina(read_reverse, read, qpos);
    // If locus is on the last position of the read and reverse, the C of the CG is not included
    // TODO: Vermeiden, wenn wir herausfinden, dass der Sequencer sich vertuen koennte.
    if c_not_included || mutation_occurred || is_invalid {
        return Some((LogProb::from(Prob(0.5)), LogProb::from(Prob(0.5))));

        return None;
    }

    Some(compute_probs_illumina(read_reverse, read, qpos))
}

/// Computes the probability of methylation/no methylation of a given position in an Illumina read. Takes mapping probability into account
///
/// # Returns
///
/// prob_alt: Probability of methylation (alternative)
/// prob_ref: Probability of no methylation (reference)
fn compute_probs_illumina(
    read_reverse: bool,
    record: &Rc<Record>,
    qpos: i32,
) -> (LogProb, LogProb) {
    let (pos_in_read, ref_base, bisulfite_base) = if !read_reverse {
        (qpos, b'C', b'T')
    } else {
        (qpos + 1, b'G', b'A')
    };

    let read_base = unsafe { record.seq().decoded_base_unchecked(pos_in_read as usize) };
    // Base quality
    let base_qual = unsafe { *record.qual().get_unchecked(pos_in_read as usize) };

    let prob_alt = prob_read_base(read_base, ref_base, base_qual);
    let no_c = if read_base != ref_base {
        read_base
    } else {
        bisulfite_base
    };
    let prob_ref = prob_read_base(read_base, no_c, base_qual);
    (prob_alt, prob_ref)
}

fn process_read_pb_np(
    read: &Rc<Record>,
    locus: &SingleLocus,
    meth_info: &Option<HashMap<usize, LogProb>>,
) -> Option<(LogProb, LogProb)> {
    let qpos = get_qpos(read, locus)?;
    // let record = &read.into_single_end_evidence()[0];
    let read_reverse = SingleLocus::read_reverse_strand(read.inner.core.flag);
    let mutation_occurred = mutation_occurred_pb_np(read_reverse, read, qpos);
    // If locus is on the last position of the read and reverse, the C of the CG is not included
    let c_not_included = qpos as usize == read.seq().len() - 1 && read_reverse;

    // If the base of the read under consideration is not a C we can't say anything about the methylation status
    if mutation_occurred || meth_info.is_none() || c_not_included {
        return None;
    }

    Some(compute_probs_pb_np(
        read_reverse,
        meth_info.as_ref().unwrap(),
        qpos,
    ))
}

/// Computes the probability of methylation/no methylation of a given position in an PacBio/ Nanopore read. Takes mapping probability into account
///
/// # Returns
///
/// prob_alt: Probability of methylation (alternative)
/// prob_ref: Probability of no methylation (reference)
fn compute_probs_pb_np(
    read_reverse: bool,
    pos_to_probs: &HashMap<usize, LogProb>,
    qpos: i32,
) -> (LogProb, LogProb) {
    let pos_in_read = qpos + if read_reverse { 1 } else { 0 };
    let prob_alt;
    let prob_ref;
    if let Some(value) = pos_to_probs.get(&(pos_in_read as usize)) {
        prob_alt = value.to_owned();
        prob_ref = LogProb::from(Prob(1_f64 - prob_alt.0.exp()));
    } else {
        prob_alt = LogProb::from(Prob(0.0));
        prob_ref = LogProb::from(Prob(1.0));
    }
    (prob_alt, prob_ref)
}

/// Computes if a mutation occured at the C of a CpG position (Not used right now but good for debugging)
///
/// # Returns
///
/// bool: True, if mutation occured
fn mutation_occurred_illumina(read_reverse: bool, record: &Rc<Record>, qpos: i32) -> bool {
    if read_reverse {
        let read_base = unsafe { record.seq().decoded_base_unchecked((qpos + 1) as usize) };
        if read_base == b'C' || read_base == b'T' {
            return true;
        }
    } else {
        let read_base = unsafe { record.seq().decoded_base_unchecked(qpos as usize) };
        if read_base == b'A' || read_base == b'G' {
            return true;
        }
    }
    return false
}

/// Computes if a given read is valid (Right now we only accept specific flags)
///
/// # Returns
///
/// bool: True, if read is valid, else false
fn read_invalid(flag: u16) -> bool {
    if flag == 0 || flag == 16 || flag == 99 || flag == 83 || flag == 147 || flag == 163 {
        return false;
    }
    true
}

/// Computes if a mutation occured at the C of a CpG position
///
/// # Returns
///
/// bool: True, if mutation occured
fn mutation_occurred_pb_np(read_reverse: bool, record: &Rc<Record>, qpos: i32) -> bool {
    if read_reverse {
        let read_base = unsafe { record.seq().decoded_base_unchecked((qpos + 1) as usize) };
        if read_base == b'C' || read_base == b'T' || read_base == b'A' {
            return true;
        }
    } else {
        let read_base = unsafe { record.seq().decoded_base_unchecked(qpos as usize) };
        if read_base == b'A' || read_base == b'G' || read_base == b'T' {
            warn!(
                "The record {:?} is not considered because a mutation occured",
                String::from_utf8_lossy(record.qname())
            );
            return true;
        }
    }
    false
}

impl Variant for Methylation {
    type Evidence = PairedEndEvidence;
    type Loci = SingleLocus;

    fn is_imprecise(&self) -> bool {
        false
    }

    /// Determine whether the evidence is suitable to assessing probabilities
    /// (i.e. overlaps the variant in the right way).
    ///
    /// # Returns
    ///
    /// The index of the loci for which this evidence is valid, `None` if invalid.
    fn is_valid_evidence(
        &self,
        evidence: &Self::Evidence,
        _: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        if match evidence {
            // TODO: Potential issues if the CpG position is at the start or end of the read.
            // This has been addressed when the C of the CpG position is at the last position of the read
            // and it's a reverse read (in which case it's not covered).
            // What still needs to be handled: when the G of a CpG position is at the first position of the read and the read is reverse.
            PairedEndEvidence::SingleEnd(read) => {
                !self.locus.overlap(read.record(), true).is_none()
            }
            PairedEndEvidence::PairedEnd { left, right } => {
                !self.locus.overlap(left.record(), true).is_none()
                    || !self.locus.overlap(right.record(), true).is_none()
                    || !self.locus.outside_overlap(left.record())
                    || !self.locus.outside_overlap(right.record())
            }
        } {
            if match self.readtype {
                // Some single PacBio reads don't have an MM:Z value and are therefore not legal
                Readtype::Illumina => true,
                Readtype::PacBio => mm_exist(&evidence.into_single_end_evidence()[0]),
                Readtype::Nanopore => mm_exist(&evidence.into_single_end_evidence()[0]),
            } {
                Some(vec![0])
            } else {
                None
            }
        } else {
            None
        }
    }

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci {
        &self.locus
    }

    fn allele_support(
        &self,
        read: &Self::Evidence,
        _alignment_properties: &AlignmentProperties,
        _alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        let qpos = match read {
            PairedEndEvidence::SingleEnd(record) => get_qpos(record.record(), &self.locus),
            PairedEndEvidence::PairedEnd { left, right } => {
                let result1 = get_qpos(left.record(), &self.locus);
                let result2 = get_qpos(right.record(), &self.locus);

                match (result1, result2) {
                    (Some(inner_qpos), _) => Some(inner_qpos),
                    (_, Some(inner_qpos)) => Some(inner_qpos),
                    _ => None,
                }
            }
        };
        if let Some(qpos) = qpos {
            let prob_alt;
            let prob_ref;
            // TODO Do something, if the next base is no G
            match read {
                PairedEndEvidence::SingleEnd(record) => match self.readtype {
                    Readtype::Illumina => {
                        (prob_alt, prob_ref) = process_read_illumina(record.record(), &self.locus)
                            .unwrap_or((LogProb(0.0), LogProb(0.0)));
                    }

                    Readtype::PacBio | Readtype::Nanopore => {
                        let meth_info = read.get_methylation_probs()[0].as_ref().unwrap_or(&None);
                        (prob_alt, prob_ref) =
                            process_read_pb_np(record.record(), &self.locus, meth_info)
                                .unwrap_or((LogProb(0.0), LogProb(0.0)));
                    }
                },

                PairedEndEvidence::PairedEnd { left, right } => {
                    match self.readtype {
                        Readtype::Illumina => {
                            // If the CpG site is only in one read included compute the probability for methylation in the read. If it is includedin both reads combine the probs.
                            let (prob_alt_left, prob_ref_left) =
                                process_read_illumina(left.record(), &self.locus)
                                    .unwrap_or((LogProb(0.0), LogProb(0.0)));
                            let (prob_alt_right, prob_ref_right) =
                                process_read_illumina(right.record(), &self.locus)
                                    .unwrap_or((LogProb(0.0), LogProb(0.0)));

                            prob_alt = LogProb(prob_alt_left.0 + prob_alt_right.0);
                            prob_ref = LogProb(prob_ref_left.0 + prob_ref_right.0);
                        }
                        // PacBio reads are normally no paired-end reads. Since we take supplementary alignments into consideration, some of the SingleEndAlignments become PairedEnd Alignments
                        // In this case we just chose the first alignment
                        Readtype::PacBio | Readtype::Nanopore => {
                            let meth_info_left =
                                read.get_methylation_probs()[0].as_ref().unwrap_or(&None);
                            let meth_info_right =
                                read.get_methylation_probs()[1].as_ref().unwrap_or(&None);

                            let (prob_alt_left, prob_ref_left) =
                                process_read_pb_np(left.record(), &self.locus, meth_info_left)
                                    .unwrap_or((LogProb(0.0), LogProb(0.0)));
                            let (prob_alt_right, prob_ref_right) =
                                process_read_pb_np(right.record(), &self.locus, meth_info_right)
                                    .unwrap_or((LogProb(0.0), LogProb(0.0)));

                            prob_alt = LogProb(prob_alt_left.0 + prob_alt_right.0);
                            prob_ref = LogProb(prob_ref_left.0 + prob_ref_right.0);
                        }
                    }
                }
            }

            let strand = if prob_ref != prob_alt {
                let record = match &read {
                    PairedEndEvidence::SingleEnd(record) => record,
                    PairedEndEvidence::PairedEnd { left, right: _ } => left,
                };

                Strand::from_record_and_pos(record.record(), qpos as usize)?
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
                    .read_position(Some(qpos as u32))
                    .third_allele_evidence(None)
                    .build()
                    .unwrap(),
            ))
        } else {
            // a read that spans an SNV might have the respective position in the
            // reference skipped (Cigar op 'N'), and the library should not choke on those reads
            // but instead needs to know NOT to add those reads (as observations) further up
            Ok(None)
        }
    }

    /// Calculate probability to sample a record length like the given one from the alt allele.
    fn prob_sample_alt(&self, _: &Self::Evidence, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

impl ToVariantRepresentation for Methylation {
    fn to_variant_representation(&self) -> model::Variant {
        model::Variant::Methylation()
    }
}

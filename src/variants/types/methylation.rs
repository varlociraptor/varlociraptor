use super::ToVariantRepresentation;
use crate::variants::evidence::realignment::Realignable;
use crate::variants::model;
use crate::{estimation::alignment_properties::AlignmentProperties, variants::sample::Readtype};

use crate::variants::evidence::bases::prob_read_base;
use crate::variants::evidence::observations::read_observation::Strand;
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, SingleEndEvidence, PairedEndEvidence, SingleLocus, Variant,
};
use rust_htslib::bam::Record;
use std::rc::Rc;
use std::collections::HashMap;

use anyhow::Result;
use bio::stats::{LogProb, Prob};
use bio_types::genome::{self, AbstractInterval, AbstractLocus};
use log::{error, warn};
use rust_htslib::bam::record::Aux;
use std::sync::Mutex;
use lazy_static::lazy_static;

// We save methylation info of the single reads for PacBio and Nanopore in order to not recompute the information for every candidate
lazy_static! {
    static ref READ_TO_METH_PROBS: Mutex<HashMap<String, HashMap<usize, f64>>> = Mutex::new(HashMap::new());
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
pub fn mm_exist(read: &SingleEndEvidence) -> bool {
    let mm_tag_exists = match (read.aux(b"Mm"), read.aux(b"MM")) {
        (Ok(_), _) | (_, Ok(_)) => true,  // True, wenn einer der Tags existiert
        _ => false,  // False, wenn keiner der Tags existiert
    };
    mm_tag_exists
}
    

/// Computes the positions of methylated Cs in PacBio and Nanopore read data
///
/// # Returns
///
/// pos_methylated_cs: Vector of positions of methylated Cs in Read
pub fn meth_pos(read: &Rc<Record>) -> Result<Vec<usize>, String> {

    let mm_tag = match (read.aux(b"Mm"), read.aux(b"MM")) {
        (Ok(tag), _) => tag,
        (_, Ok(tag)) => tag,
        _ => {
            warn!("MM value not found on pos {:?}", read.inner.core.pos);
            return Ok(Vec::new());

        }
    };
    
    // let mm_tag = read.aux(b"Mm").map_err(|e| e.to_string())?;
    if let Aux::String(tag_value) = mm_tag {
        let mut mm = tag_value.to_owned();
        
        if !mm.is_empty() {
            let read_reverse = read_reverse_strand(read.inner.core.flag);
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
            if let Some(methylated_part) = mm.strip_prefix("C+m,") {
                let mut meth_pos = 0;
                
                let mut pos_methylated_cs: Vec<usize> = methylated_part
                    .split(',')
                    .filter_map(|position_str| {
                        position_str.parse::<usize>().ok().map(|position| {
                            meth_pos += position + 1;
                            meth_pos
                        }) // methylated Cs
                    })
                    .map(|pos| pos_read_base[pos - 1])// Chose only the positions of methylated Cs out of all Cs
                    .collect();

                if read_reverse {
                    pos_methylated_cs.reverse();
                } 
                return Ok(pos_methylated_cs);
            }
        }
        // No methylation info in this read
        else {
            return Ok(Vec::new());
        }
    } else {
        warn!("MM tag in bam file is not valid on pos {:?}", read.inner.core.pos);
        return Ok(Vec::new());
    }
    warn!("Error while obtaining MM:Z tag on pos {:?}", read.inner.core.pos);
    Ok(Vec::new())
}


/// Computes the probability of methylation in PacBio and Nanopore read data
///
/// # Returns
///
/// ml: Vector of methylation probabilities
pub fn meth_probs(read: &Rc<Record>) -> Result<Vec<LogProb>, String> {
    let ml_tag = match (read.aux(b"Ml"), read.aux(b"ML")) {
        (Ok(tag), _) => tag,
        (_, Ok(tag)) => tag,
        _ => {
            warn!("ML value not found on pos {:?}", read.inner.core.pos);
            return Ok(Vec::new());

        }
    };
    // let ml_tag = read.aux(b"Ml").map_err(|e| e.to_string())?;
    let read_reverse = read_reverse_strand(read.inner.core.flag);
    if let Aux::ArrayU8(tag_value) = ml_tag {
        let mut ml: Vec<LogProb> = tag_value.iter().map(|val| LogProb::from(Prob((f64::from(val) + 0.5) / 256.0))).collect();
        if read_reverse{
            ml.reverse();
        }
        Ok(ml)
    } else {
        warn!("MM tag in bam file is not valid on pos {:?}", read.inner.core.pos);
        Ok(Vec::new())
    }
}

/// Computes the probability of methylation/no methylation of a given position in an PacBio/ Nanopore read. Takes mapping probability into account
///
/// # Returns
///
/// prob_alt: Probability of methylation (alternative)
/// prob_ref: Probability of no methylation (reference)
fn compute_probs_pb_np(pos_in_read: i32, pos_to_probs: &HashMap<usize, LogProb>) -> (LogProb, LogProb) {
    let prob_alt;
    let prob_ref;
    if let Some(value) = pos_to_probs.get(&(pos_in_read as usize)) {
        prob_alt = value.to_owned();
        prob_ref = LogProb::from(Prob(1 as f64 - prob_alt.0.exp()));
    } else {
        prob_alt = LogProb::from(Prob(0.01));
        prob_ref = LogProb::from(Prob(0.99));
    }
    (prob_alt, prob_ref)
}


/// Position of CpG site in the read
///
/// # Returns
///
/// Option<position>: Position of CpG site if it is included in read
/// None else
fn get_qpos(read: &Rc<Record> , locus: &SingleLocus) -> Option<i32> {
    read
        .cigar_cached()
        .unwrap()
        .read_pos(locus.range().start as u32, false, false).unwrap().map(|qpos| qpos as i32)
}

/// Computes the probability of methylation/no methylation of a given position in an Illumina read. Takes mapping probability into account
///
/// # Returns
///
/// prob_alt: Probability of methylation (alternative)
/// prob_ref: Probability of no methylation (reference)
fn compute_probs_illumina(reverse_read: bool, record:  &Rc<Record>, qpos: i32) -> (LogProb, LogProb){
    let (pos_in_read, ref_base, bisulfite_base) = if !reverse_read {
        (qpos, b'C', b'T')
    } else {
        (qpos + 1, b'G', b'A')
    };
          
    let read_base = unsafe { record.seq().decoded_base_unchecked(pos_in_read as usize) };
    let base_qual = unsafe { *record.qual().get_unchecked(pos_in_read as usize) };

    let prob_alt = prob_read_base(read_base, ref_base, base_qual);
    let no_c = if read_base != ref_base { read_base } else { bisulfite_base };
    let prob_ref = prob_read_base(read_base, no_c, base_qual);
    (prob_alt, prob_ref)
}


fn process_read(read: &Rc<Record>, locus: &SingleLocus) -> Option<(LogProb, LogProb)> {
    let qpos = get_qpos(read, locus)?;
    let reverse_read = read_reverse_strand(read.inner.core.flag);
    let is_invalid = read_invalid(read.inner.core.flag);
    let mutation_occurred = mutation_occurred(reverse_read, read, qpos);

    // If locus is on the last position of the read and reverse, the C of the CG is not included
    if (qpos as usize == read.seq().len() - 1 && reverse_read) || is_invalid || mutation_occurred {
        return None;
    }

    Some(compute_probs_illumina(reverse_read, read, qpos))
}

/// Finds out whether the given string is a forward or reverse string.
///
/// # Returns
///
/// True if read given read is a reverse read, false if it is a forward read
pub fn read_reverse_strand(flag:u16) -> bool {
    let read_paired = 0b1;
    let read_mapped_porper_pair = 0b01;
    let read_reverse = 0b10000;
    let mate_reverse = 0b100000;
    let first_in_pair = 0b1000000;
    let second_in_pair = 0b10000000;
    if  (flag & read_paired != 0 && flag & read_mapped_porper_pair != 0 && flag & read_reverse != 0 && flag & first_in_pair != 0) ||
        (flag & read_paired != 0 && flag & read_mapped_porper_pair != 0 && flag & mate_reverse != 0 && flag & second_in_pair != 0) ||
        (flag & read_reverse != 0 && (flag & read_paired == 0 || flag & read_mapped_porper_pair == 0)) {
            return true;
    }
    false
}

/// Computes if a given read is valid (Right now we only accept specific flags)
///
/// # Returns
///
/// bool: True, if read is valid, else false
fn read_invalid(flag:u16) -> bool {
    if flag == 0 || flag == 16 || flag == 99 || flag == 83 || flag == 147 || flag == 163 {
        return false
    }
    true
}

/// Computes if a mutation occured at the C of a CpG position
///
/// # Returns
///
/// bool: True, if mutation occured
fn mutation_occurred(reverse_read: bool, record:  &Rc<Record>, qpos: i32) -> bool {
    if reverse_read {
        let read_base = unsafe { record.seq().decoded_base_unchecked((qpos + 1) as usize) };
        if read_base == b'C' || read_base == b'T' {
            return  true
        }
    }
    else {
        let read_base = unsafe { record.seq().decoded_base_unchecked(qpos as usize) };
        if read_base == b'A' || read_base == b'G' {
            return  true
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
            // TODO maybe problems, if CpG position starts just at the beginning/end of the read (Ich habe mich eigentlich darum gekuemmert, wenn das C der CpG Position die letzte Stelle des Reads is und es ein reverse Read ist (dann ist es nicht abgesdeckt.) Was fehlt: Wenn die erste Position das G einer CpG Stelle ist und der Read reverse ist.)
            PairedEndEvidence::SingleEnd(read) => !self.locus.overlap(read.record(), true).is_none(),
            PairedEndEvidence::PairedEnd { left, right } => {
                !self.locus.overlap(left.record(), true).is_none()
                    || !self.locus.overlap(right.record(), true).is_none() 
                    || !self.locus.outside_overlap(left.record())
                    || !self.locus.outside_overlap(right.record())
            }
        } {
            if match self.readtype {
                    // Some single PacBio reads don't have an MM:Z value and are therefore not legal
                    Readtype::Illumina => {true},
                    Readtype::PacBio => {mm_exist(&evidence.into_single_end_evidence()[0])},
                    Readtype::Nanopore => {mm_exist(&evidence.into_single_end_evidence()[0])},

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
            PairedEndEvidence::SingleEnd(record) => {
                get_qpos(record.record(), &self.locus)
            }
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
                PairedEndEvidence::SingleEnd(record) => {
                    match self.readtype {

                        Readtype::Illumina => {
                            let reverse_read = read_reverse_strand(record.record().inner.core.flag);
                            if mutation_occurred(reverse_read, record.record(), qpos) || read_invalid(record.record().inner.core.flag) {
                                return Ok(None)
                            }
                            (prob_alt, prob_ref) = compute_probs_illumina(reverse_read, record.record(), qpos);
                        }

                        Readtype::PacBio | Readtype::Nanopore=> {
                            let record = &read.into_single_end_evidence()[0];  
                            let read_reverse = read_reverse_strand(record.inner.core.flag);
                            let pos_in_read = qpos + if read_reverse { 1 } else { 0 };
                            let read_base = unsafe { record.seq().decoded_base_unchecked((pos_in_read) as usize) };  
                            // If the base of the read under consideration is not a C we can't say anything about the methylation status  
                            if !((read_base == b'C' && !read_reverse) || (read_base == b'G' && read_reverse)) {    
                                return Ok(None)
                            }
                            (prob_alt, prob_ref) = compute_probs_pb_np(pos_in_read, read.get_methylation_probs()[0].as_ref().unwrap());  
                        }
                    } 
                }      
                                 
                PairedEndEvidence::PairedEnd { left, right } => {
                    match self.readtype {

                        Readtype::Illumina => {
                            // If the CpG site is only in one read included compute the probability for methylation in the read. If it is includedin both reads combine the probs.
                            let (prob_alt_left, prob_ref_left) = process_read(left.record(), &self.locus).unwrap_or((LogProb(0.0), LogProb(0.0)));
                            let (prob_alt_right, prob_ref_right) = process_read(right.record(), &self.locus).unwrap_or((LogProb(0.0), LogProb(0.0)));
                            
                            prob_alt = LogProb(prob_alt_left.0 + prob_alt_right.0);
                            prob_ref = LogProb(prob_ref_left.0 + prob_ref_right.0);
                        }
                        // PacBio reads are normally no paired-end reads. Since we take supplementary alignments into consideration, some of the SingleEndAlignments become PairedEnd Alignments
                        // In this case we just chose the first alignment
                        Readtype::PacBio | Readtype::Nanopore => {
                            let record = &read.into_single_end_evidence()[0];  
                            let read_reverse = read_reverse_strand(record.inner.core.flag);
                            let pos_in_read = qpos + if read_reverse { 1 } else { 0 };
                            let read_base = unsafe { record.seq().decoded_base_unchecked((pos_in_read) as usize) };                              
                            // If the base of the read under consideration is not a C we can't say anything about the methylation status  
                            if !((read_base == b'C' && !read_reverse) || (read_base == b'G' && read_reverse)) {    
                                return Ok(None)
                            }
                            (prob_alt, prob_ref) = compute_probs_pb_np(pos_in_read, read.get_methylation_probs()[0].as_ref().expect("No meth value"));  
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

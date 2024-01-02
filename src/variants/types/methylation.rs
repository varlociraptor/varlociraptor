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
    static ref READ_TO_METH_PROBS: Mutex<HashMap<i64, HashMap<usize, f64>>> = Mutex::new(HashMap::new());
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

/// With MM tags it is important to know whether it is a forward or reverse strand. With forward, you go from front to back and count the Cs. With reverse, you go from back to front in the sequence, in the MM tag from front to back and count the Gs
fn meth_pos(read: &SingleEndEvidence) -> Result<Vec<usize>, String> {
    let mm_tag = read.aux(b"Mm").map_err(|e| e.to_string())?;
    let read_reverse = read_reverse_strand(read.inner.core.flag);
    if let Aux::String(tag_value) = mm_tag {
        let mut mm = tag_value.to_owned();
        if !mm.is_empty() {
            // Compute the positions of all Cs (or Gs for reverse strand) in the read
            let read_seq = String::from_utf8_lossy(&read.seq().as_bytes()).to_string();
            let pos_read_base: Vec<usize>;
            if !read_reverse{
                pos_read_base = read_seq
                    .char_indices()
                    .filter(|(_, c)| *c == 'C')
                    .map(|(index, _)| index)
                    .collect();
            }
            else {
                pos_read_base = read_seq
                .char_indices()
                .filter(|(_, c)| *c == 'G')
                .map(|(index, _)| index)
                .rev()
                .collect();
            }
            // Compute which Cs are methylated
            mm.pop();
            if let Some(methylated_part) = mm.strip_prefix("C+m,") {
                let mut meth_pos = 0;
                let mut methylated_cs: Vec<usize> = methylated_part
                    .split(',')
                    .filter_map(|position_str| {
                        position_str.parse::<usize>().ok().map(|position| {
                            meth_pos += position + 1;
                            meth_pos
                        })
                    })
                    .collect();
                // If last C is not methylated, there has been added one C to much
                if methylated_cs[methylated_cs.len() - 1] > pos_read_base.len() {
                    methylated_cs.pop();
                }
                // Chose only the methylated Cs out of all Cs
                let mut pos_methylated_cs: Vec<usize> = 
                    methylated_cs.iter().map(|&pos| pos_read_base[pos - 1]).collect();
            
                if read_reverse {
                    pos_methylated_cs.reverse();
                } 
                return Ok(pos_methylated_cs);
            }
        }
    } else {
        error!("MM tag in bam file is not valid");
    }
    Err("Error while obtaining MM:Z tag".to_string())
}


/// Computes the probability of methylation in PacBio read data
///
/// # Returns
///
/// prob_alt: Probability of methylation (alternative)
/// prob_ref: Probaability of no methylation (reference)
fn meth_probs(read: &SingleEndEvidence) -> Result<Vec<f64>, String> {
    let ml_tag = read.aux(b"Ml").map_err(|e| e.to_string())?;
    let read_reverse = read_reverse_strand(read.inner.core.flag);
    if let Aux::ArrayU8(tag_value) = ml_tag {
        let mut ml: Vec<f64> = tag_value.iter().map(|val| f64::from(val) / 255.0).collect();
        if read_reverse{
            ml.reverse();
        }
        return Ok(ml);

    } else {
        error!("Tag is not of type String");
    }
    Err("Error while obtaining ML:B tag".to_string())
}



/// Position of CpG site in the read
///
/// # Returns
///
/// Option<position> if CpG site included in read
/// None else
fn get_qpos(read: &Rc<Record> , locus: &SingleLocus) -> Option<i32> {
    if let Some(qpos) = read
        .cigar_cached()
        .unwrap()
        .read_pos(locus.range().start as u32, false, false).unwrap()
    {
        Some(qpos as i32)
    } else {
        None
    }
}

/// Computes the probability of methylation/no methylation of a given position in a read. Takes mapping probability into account
///
/// # Returns
///
/// prob_alt: Probability of methylation (alternative)
/// prob_ref: Probaability of no methylation (reference)
fn compute_probs(reverse_read: bool, record:  &Rc<Record>, qpos: i32) -> (LogProb, LogProb){
    let prob_alt;
    let prob_ref;
    if !reverse_read {          
        let read_base = unsafe { record.seq().decoded_base_unchecked(qpos as usize) };
        let base_qual = unsafe { *record.qual().get_unchecked(qpos as usize) };
        // Prob_read_base: Wkeit, dass die gegebene Readbase tatsachlich der 2. base entspricht (Also dass es eigtl die 2. Base ist)
        prob_alt = prob_read_base(read_base, b'C', base_qual);
        let no_c = if read_base != b'C' { read_base } else { b'T' };
        prob_ref = prob_read_base(read_base, no_c, base_qual);
    }
    else {
        // In the reverse String we have to check on the position qpos + 1, because CG ist in the reverse string GC.
        // warn!("Qpos: {:?}, Sequence: {:?}", qpos, record.seq().len());

        let read_base = unsafe { record.seq().decoded_base_unchecked((qpos + 1) as usize) };      
        let base_qual = unsafe { *record.qual().get_unchecked((qpos + 1) as usize) };
        // After aligning every base is flipped so we want the probabilitz for G and not for C (now GC turns into CG)
        prob_alt = prob_read_base(read_base, b'G', base_qual);
        let no_g = if read_base != b'G' { read_base } else { b'A' };
        prob_ref = prob_read_base(read_base, no_g, base_qual);
    }
    (prob_alt, prob_ref)
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
    if (flag & read_paired) != 0 && (flag & read_mapped_porper_pair) != 0 {
        if (flag & read_reverse) != 0 && (flag & first_in_pair) != 0 {
            return true
        }
        else if (flag & mate_reverse) != 0 && (flag & second_in_pair) != 0 {
            return true
        }
    }
    else {
        if (flag & read_reverse) != 0 {
            return true
        }
    }
    false
    // read.inner.core.flag == 163 || read.inner.core.flag == 83 || read.inner.core.flag == 16
}

fn read_invalid(flag:u16) -> bool {
    // let secondary_alignment = 0b100000000;
    // let qc_failed = 0b1000000000;
    // let duplicate = 0b10000000000;
    // let supplemental = 0b100000000000;
    // if (flag & secondary_alignment) != 0 || (flag & qc_failed) != 0 
    // || (flag & duplicate) != 0 || (flag & supplemental) != 0 {
    //     return true
    // }
    // false    
    if flag == 0 || flag == 16 || flag == 99 || flag == 83 || flag == 147 || flag == 163 {
        return false
    }
    return true
}

fn mutation_occurred(reverse_read: bool, record:  &Rc<Record>, qpos: i32) -> bool {
    if reverse_read {
        let read_base = unsafe { record.seq().decoded_base_unchecked((qpos + 1) as usize) };
        if read_base == 67 || read_base == 84 {
            return  true
        }
    }
    else {
        let read_base = unsafe { record.seq().decoded_base_unchecked(qpos as usize) };
        if read_base == 65 || read_base == 71 {
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
            PairedEndEvidence::SingleEnd(read) => !self.locus.overlap(read, true).is_none(),
            PairedEndEvidence::PairedEnd { left, right } => {
                !self.locus.overlap(left, true).is_none()
                    || !self.locus.overlap(right, true).is_none() 
                    || !self.locus.outside_overlap(left)
                    || !self.locus.outside_overlap(right)
            }
        } {
            if match self.readtype {
                    // Some single PacBio reads don't have an MM:Z value and are therefore not legal
                    Readtype::Illumina => {true},
                    Readtype::PacBio => {!meth_pos(&evidence.into_single_end_evidence()[0]).is_err()},
                    Readtype::Nanopore => {!meth_pos(&evidence.into_single_end_evidence()[0]).is_err()},

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
        let qpos;
        qpos = match read {
            PairedEndEvidence::SingleEnd(record) => {
                    get_qpos(record, &self.locus)
            }
            PairedEndEvidence::PairedEnd { left, right } => {
                let result1 = get_qpos(left, &self.locus);
                let result2 = get_qpos(right, &self.locus);
        
                match (result1, result2) {
                    (Some(inner_qpos), _) => Some(inner_qpos),
                    (_, Some(inner_qpos)) => Some(inner_qpos),
                    _ => None,
                }
            }
        };
        if let Some(qpos) = qpos {
            let mut prob_alt = LogProb::from(Prob(0.0));
            let mut prob_ref = LogProb::from(Prob(0.0));
            // TODO Do something, if the next base is no G

            
            match read {
                PairedEndEvidence::SingleEnd(record) => {
                    match self.readtype {
                        Readtype::Illumina => {
                            // return Ok(None);
                            // TODO in paired end data are some single end reads?
                            if let Some(qpos) = get_qpos(record, &self.locus) {
                                let reverse_read = read_reverse_strand(record.inner.core.flag);
                                if mutation_occurred(reverse_read, record, qpos) || read_invalid(record.inner.core.flag) {
                                    return Ok(None)
                                }
                                (prob_alt, prob_ref) = compute_probs(reverse_read, record, qpos);
                                if self.locus.interval.range().start == 3 {
                                    warn!("Single: {:?}, {:?}", record.inner.core.flag,  record.inner.core.pos + 1);
                                    warn!("Prob_alt: {:?}", prob_alt);
                                    warn!("Prob_ref: {:?}", prob_ref);
                                }
                            }  
                        }
                        Readtype::PacBio | Readtype::Nanopore=> {
                            let record = read.into_single_end_evidence();  
                            let mut pos_in_read = qpos;
                            if read_reverse_strand(record[0].inner.core.flag) {
                                pos_in_read += 1;
                            }  
                            let mut data = READ_TO_METH_PROBS.lock().unwrap(); // 
                            // Wenn dieser Read fuer einen anderen Kandidaten bereits betrachtet wurde haben wir die MM und ML Informationen bereitsgespeichert
                            let pos_to_probs: HashMap<usize, f64>;
                            if let Some(pos_to_probs_found) = data.get(&record[0].inner.core.pos) {
                                pos_to_probs = pos_to_probs_found.clone();
                            }
                            else {
                                println!("{:?}", record[0].inner.core.pos);
                                let meth_pos = meth_pos(&record[0]).unwrap();
                                let meth_probs = meth_probs(&record[0]).unwrap();
                                pos_to_probs = meth_pos.into_iter().zip(meth_probs.into_iter()).collect();
                                data.insert(record[0].inner.core.pos, pos_to_probs.clone());                     
                            }
                         
                                                
                            if let Some(value) = pos_to_probs.get(&(pos_in_read as usize)) {
                                prob_alt = LogProb::from(Prob(*value as f64));
                                prob_ref = LogProb::from(Prob(1 as f64 - *value as f64));
                            
                            
                            } else {
                                // TODO What should I do if there is no prob given
                                prob_alt = LogProb::from(Prob(0.0));
                                prob_ref = LogProb::from(Prob(1.0));
                                // warn!("No probability given for unmethylated Cs!");
                            }
                        }
                    } 
                }                       
                PairedEndEvidence::PairedEnd { left, right } => {
                    match self.readtype {
                        Readtype::Illumina => {
                            // If the CpG site is only in one read included compute the probability for methylation in the read. If it is includedin both reads combine the probs.
                            let qpos_left = get_qpos(left, &self.locus);
                            let qpos_right = get_qpos(right, &self.locus);
                            let mut prob_alt_left = LogProb(0.0);
                            let mut prob_ref_left = LogProb(0.0);
                            let mut prob_alt_right = LogProb(0.0);
                            let mut prob_ref_right = LogProb(0.0);
                            let left_invalid = read_invalid(left.inner.core.flag);
                            let right_invalid = read_invalid(right.inner.core.flag);
                            if left_invalid && right_invalid {
                                return Ok(None);
                            }
                            if let Some(qpos_left) = qpos_left {
                                let reverse_read = read_reverse_strand(left.inner.core.flag);
                                // Don't lookat read if it is a mutation and not methylation
                                // if mutation_occurred(reverse_read, left, qpos_left) {
                                //     return Ok(None)
                                // }
                                // If locus is on last position of the read and reverse, the C of the CG is not included
                                if (qpos_left as usize == left.seq().len() - 1 && reverse_read) || left_invalid{
                                    return Ok(None)
                                }
                                (prob_alt_left, prob_ref_left) = compute_probs(reverse_read, left, qpos_left);
                            }
                            if let Some(qpos_right) = qpos_right {
                                let reverse_read = read_reverse_strand(right.inner.core.flag);
                                // Don't lookat read if it is a mutation and not methylation
                                // if mutation_occurred(reverse_read, right, qpos_right) {
                                //     return Ok(None)
                                // }
                                if (qpos_right as usize == right.seq().len() - 1 && reverse_read) || right_invalid{
                                    return Ok(None)
                                }
                                (prob_alt_right, prob_ref_right) = compute_probs(reverse_read, right, qpos_right);
                            }                                 
                            prob_alt = LogProb(prob_alt_left.0 + prob_alt_right.0);
                            prob_ref = LogProb(prob_ref_left.0 + prob_ref_right.0);
                            if self.locus.interval.range().start ==  3 {
                                warn!("Positions: {:?}, {:?}", left.inner.core.pos + 1, right.inner.core.pos + 1);

                                warn!("Left: {:?}, {:?}", left.inner.core.flag, qpos_left);
                                warn!("Right: {:?}, {:?}", right.inner.core.flag, qpos_right);
                                warn!("Prob_alt: {:?}", prob_alt);
                                warn!("Prob_ref: {:?}", prob_ref);

                            }
                        }
                        // PacBio reads are normally no paired-end reads. Since we take supplementary alignments into consideration, some of the SingleEndAlignments become PairedEnd Alignments
                        Readtype::PacBio | Readtype::Nanopore => {
                            let mut record = read.into_single_end_evidence();  
                            // Chose non-supplementary alignment
                            if left.inner.core.flag < 2000 {

                                record.remove(1);
                            }
                            else {
                                record.remove(0);

                            }



                            let record = read.into_single_end_evidence();  
                            let mut pos_in_read = qpos;
                            if read_reverse_strand(record[0].inner.core.flag) {
                                pos_in_read += 1;
                            }  
                            let mut data = READ_TO_METH_PROBS.lock().unwrap();
                            // If this read has already been viewed for another candidate, we have already saved the MM and ML information
                            let pos_to_probs: HashMap<usize, f64>;
                            if let Some(pos_to_probs_found) = data.get(&record[0].inner.core.pos) {
                                pos_to_probs = pos_to_probs_found.clone();
                            }
                            else {
                                println!("{:?}", record[0].inner.core.pos);

                                let meth_pos = meth_pos(&record[0]).unwrap();
                                let meth_probs = meth_probs(&record[0]).unwrap();
                                pos_to_probs = meth_pos.into_iter().zip(meth_probs.into_iter()).collect();
                                data.insert(record[0].inner.core.pos, pos_to_probs.clone());                     
                            }
                         
                                                
                            if let Some(value) = pos_to_probs.get(&(pos_in_read as usize)) {
                                prob_alt = LogProb::from(Prob(*value as f64));
                                prob_ref = LogProb::from(Prob(1 as f64 - *value as f64));
                            
                            
                            } else {
                                // TODO What should I do if there is no prob given
                                prob_alt = LogProb::from(Prob(0.0));
                                prob_ref = LogProb::from(Prob(1.0));
                                warn!("No probability given for unmethylated Cs!");
                            }
                        } 
                    }

                
                }
            }
            
          
            let strand = if prob_ref != prob_alt {
                    let record = match &read {
                        PairedEndEvidence::SingleEnd(record) => record,
                        PairedEndEvidence::PairedEnd { left, right: _ } => left,
                    };
                    
                    Strand::from_record_and_pos(&record, qpos as usize)?
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

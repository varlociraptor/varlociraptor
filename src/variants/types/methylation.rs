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

fn meth_pos(read: &SingleEndEvidence) -> Result<Vec<usize>, String> {
    let mm_tag = read.aux(b"MM").map_err(|e| e.to_string())?;
    if let Aux::String(tag_value) = mm_tag {
        let mut mm = tag_value.to_owned();
        if !mm.is_empty() {
            // Compute the positions of all Cs in the Read

            let read_seq = String::from_utf8_lossy(&read.seq().as_bytes()).to_string();
            let pos_cs: Vec<usize> = read_seq
                .char_indices()
                .filter(|(_, c)| *c == 'C')
                .map(|(index, _)| index)
                .collect();
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
                if methylated_cs[methylated_cs.len() - 1] > pos_cs.len() {
                    methylated_cs.pop();
                }
                // Chose only the methylated Cs out of all Cs
                let pos_methylated_cs: Vec<usize> =
                    methylated_cs.iter().map(|&pos| pos_cs[pos - 1]).collect();
                return Ok(pos_methylated_cs);
            }
        }
    } else {
        error!("Tag is not of type String");
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
    let ml_tag = read.aux(b"ML").map_err(|e| e.to_string())?;
    if let Aux::ArrayU8(tag_value) = ml_tag {
        let ml: Vec<f64> = tag_value.iter().map(|val| f64::from(val) / 255.0).collect();
        return Ok(ml);
    } else {
        error!("Tag is not of type String");
    }
    Err("Error while obtaining ML:B tag".to_string())
}

fn qpos_in_read(qpos: i32, read: &Rc<Record>) -> bool {
    return read.inner.core.pos <= qpos as i64 && qpos <= read.inner.core.pos as i32 + read.inner.core.l_qseq;
}

/// Postion of CpG site in the read
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
fn read_reverse_strand(flag:u16, paired: bool) -> bool {
    let read_reverse = 0b10000;
    let mate_reverse = 0b100000;
    let first_in_pair = 0b1000000;
    let second_in_pair = 0b10000000;
    if paired{
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
            // TODO maybe problems, if CpG position starts/ends just at the beginning/end of the read
            PairedEndEvidence::SingleEnd(read) => !self.locus.overlap(read, true).is_none(),
            PairedEndEvidence::PairedEnd { left, right } => {
                !self.locus.overlap(left, true).is_none()
                    || !self.locus.overlap(right, true).is_none()
            }
        } {
            Some(vec![0])
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
            if self.locus.interval.range().start == 8447241 {
                // warn!("Readposition: {:?} \n", read.PairedEnd.inner.core.pos);
                warn!("##############################");
                warn!("QPos: {:?}", qpos);
            }
            match self.readtype {
                Readtype::Illumina => {
                    match read {
                        PairedEndEvidence::SingleEnd(record) => {
                            // TODO in paired end data are some single end reads?
                            if let Some(qpos) = get_qpos(record, &self.locus) {
                                let reverse_read = read_reverse_strand(record.inner.core.flag, false);
                                (prob_alt, prob_ref) = compute_probs(reverse_read, record, qpos);
                            }  
                        }                        
                        PairedEndEvidence::PairedEnd { left, right } => {
                            // If the CpG site is only in one read included compute the probability for methylation in the read. If it is includedin both reads combine the probs.
                            let qpos_left = get_qpos(left, &self.locus);
                            let qpos_right = get_qpos(right, &self.locus);
                            let mut prob_alt_left = LogProb(0.0);
                            let mut prob_ref_left = LogProb(0.0);
                            let mut prob_alt_right = LogProb(0.0);
                            let mut prob_ref_right = LogProb(0.0);
                            if let Some(qpos_left) = qpos_left {
                                let reverse_read = read_reverse_strand(left.inner.core.flag, true);
                                (prob_alt_left, prob_ref_left) = compute_probs(reverse_read, left, qpos_left);
                            }
                            if let Some(qpos_right) = qpos_right {
                                let reverse_read = read_reverse_strand(right.inner.core.flag, true);
                                (prob_alt_right, prob_ref_right) = compute_probs(reverse_read, right, qpos_right);
                            }                                 
                            prob_alt = LogProb(prob_alt_left.0 + prob_alt_right.0);
                            prob_ref = LogProb(prob_ref_left.0 + prob_ref_right.0);
                        }
                    }
                }

                Readtype::PacBio => {
                    // prob_alt = LogProb::from(Prob(0.0));
                    // prob_ref = LogProb::from(Prob(1.0));
                    // warn!("PacBio");
                    // let record = read.into_single_end_evidence();  
                    //  // // Get methylation info from MM and ML TAG.
                    // let meth_pos = meth_pos(record).unwrap();
                    // let meth_probs = meth_probs(read).unwrap();
                    // let pos_to_probs: HashMap<usize, f64> =
                    //     meth_pos.into_iter().zip(meth_probs.into_iter()).collect();
                    // if let Some(value) = pos_to_probs.get(&(qpos as usize)) {
                    //     prob_alt = LogProb::from(Prob(*value as f64));
                    //     prob_ref = LogProb::from(Prob(1 as f64 - *value as f64));
                    // } else {
                    //     // TODO What should I do if there is no prob given
                    //     prob_alt = LogProb::from(Prob(0.0));
                    //     prob_ref = LogProb::from(Prob(1.0));
                    //     warn!("No probability given for unmethylated Cs!");
                    // }
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

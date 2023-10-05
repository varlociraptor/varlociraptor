// cargo run -- preprocess variants ~/Documents/Promotion/varlociraptor-methylation-evaluation/resources/example-genome.fasta --candidates ~/Documents/Promotion/varlociraptor-methylation-evaluation/resources/example-candidates.bcf --bam ~/Documents/Promotion/varlociraptor-methylation-evaluation/resources/example-reads.bam > ~/Documents/Promotion/varlociraptor-methylation-evaluation/resources/observations.bcf

use super::ToVariantRepresentation;
use crate::variants::evidence::realignment::Realignable;
use crate::variants::model;
use crate::{estimation::alignment_properties::AlignmentProperties, variants::sample::Readtype};

use crate::variants::evidence::bases::prob_read_base;
use crate::variants::evidence::observations::read_observation::Strand;
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, Overlap, MultiLocus, SingleEndEvidence, PairedEndEvidence, SingleLocus, Variant,
};
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
    
    // #################################################################### - Old SingleEnd - ###################################################
    // fn is_valid_evidence(
    //     &self,
    //     evidence: &SingleEndEvidence,
    //     _: &AlignmentProperties,
    // ) -> Option<Vec<usize>> {
    //     if let Overlap::Enclosing = self.locus.overlap(evidence, false) {
    //         Some(vec![0])
    //     // If the forward read end with a C of a CpG site or the reverse reag starts with a C in a CpG site include the read
    //     } else if self.locus.outside_overlap(evidence) {
    //         Some(vec![0])
    //     } else {
    //         None
    //     }
    // }
    // #################################################################### - Old SingleEnd - ###################################################



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
        // qpos: Position of the C under consideration in the read 
        // let reverse_read =  (read.inner.core.flag & 0x10) != 0; // If the Flag Contains 16 (in hex 0x10), the read is a reverse read
        // let reverse_read = (read.inner.core.flag & 0x10) != 0;
        // let reverse_read = read.inner.core.flag == 163 || read.inner.core.flag == 83 || read.inner.core.flag == 16;
        let reverse_read = true;
        let qpos = Some(10);
        // let qpos = if let Some(inner_qpos) = read
        //     .cigar_cached()
        //     .unwrap()
        //     .read_pos(self.locus.range().start as u32, false, false)?
        // {
        //     Some(inner_qpos as i32)
        // } else if self.locus.outside_overlap(read) {
        //     Some(if reverse_read { -1 } else { (read.cigar_len() + 1) as i32 })
        // } else {
        //     None
        // };
            
        if let Some(qpos) = qpos {
            let prob_alt;
            let prob_ref;
            // TODO Do something, if the next base is no G
            match self.readtype {
                Readtype::Illumina => {
                    match read {
                        PairedEndEvidence::SingleEnd(record) => {
                            if !reverse_read {
                                // let read_base = unsafe { read.seq().decoded_base_unchecked(qpos as usize) };
                                // let base_qual = unsafe { *read.qual().get_unchecked(qpos as usize) };
                                // // Prob_read_base: Wkeit, dass die gegebene Readbase tatsachlich der 2. base entspricht (Also dass es eigtl die 2. Base ist)
                                // prob_alt = prob_read_base(read_base, b'C', base_qual);
                                // let no_c = if read_base != b'C' { read_base } else { b'T' };
                                // prob_ref = prob_read_base(read_base, no_c, base_qual);
                                println!("Single und falschrum")
                            }
                            else {
                                // let read_base = unsafe { read.seq().decoded_base_unchecked((qpos + 1) as usize) };
                                // let base_qual = unsafe { *read.qual().get_unchecked((qpos + 1) as usize) };
                                // prob_alt = prob_read_base(read_base, b'G', base_qual);
                                // let no_g = if read_base != b'G' { read_base } else { b'A' };
                                // prob_ref = prob_read_base(read_base, no_g, base_qual);
                                println!("Single und richtigrum")

                            }
                        }
                        PairedEndEvidence::PairedEnd { left, right } => {
                            println!("Okay, hier sind tatsaechlich paired end Reads!");
                        }
                    }
                }
                Readtype::PacBio => {
                    // // Get methylation info from MM and ML TAG.
                    // let meth_pos = meth_pos(read).unwrap();
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
            // TODO: Implement strand
            // let strand = if prob_ref != prob_alt {
            //     Strand::from_record_and_pos(read, qpos as usize)?
            // } else {
            //     // METHOD: if record is not informative, we don't want to
            //     // retain its information (e.g. strand).
            //     Strand::no_strand_info()
            // };
            let strand = Strand::no_strand_info();
            Ok(Some(
                AlleleSupportBuilder::default()
                    .prob_ref_allele(prob_ref)
                    .prob_alt_allele(prob_alt)
                    .strand(strand)
                    .read_position(Some(qpos as u32))
                    // TODO: Implement third allele
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
 // #################################################################### - Old SingleEnd - ###################################################
    // fn allele_support(
    //     &self,
    //     read: &SingleEndEvidence,
    //     _alignment_properties: &AlignmentProperties,
    //     _alt_variants: &[Box<dyn Realignable>],
    // ) -> Result<Option<AlleleSupport>> {
    //     // qpos: Position of the C under consideration in the read 
    //     // let reverse_read =  (read.inner.core.flag & 0x10) != 0; // If the Flag Contains 16 (in hex 0x10), the read is a reverse read
    //     // let reverse_read = (read.inner.core.flag & 0x10) != 0;
    //     let reverse_read = read.inner.core.flag == 163 || read.inner.core.flag == 83 || read.inner.core.flag == 16;
    //     let qpos = if let Some(inner_qpos) = read
    //         .cigar_cached()
    //         .unwrap()
    //         .read_pos(self.locus.range().start as u32, false, false)?
    //     {
    //         Some(inner_qpos as i32)
    //     } else if self.locus.outside_overlap(read) {
    //         Some(if reverse_read { -1 } else { (read.cigar_len() + 1) as i32 })
    //     } else {
    //         None
    //     };
            
    //     if let Some(qpos) = qpos {
    //         let prob_alt;
    //         let prob_ref;
    //         // TODO Do something, if the next base is no G
    //         match self.readtype {
    //             Readtype::Illumina => {
    //                 if !reverse_read {
    //                     let read_base = unsafe { read.seq().decoded_base_unchecked(qpos as usize) };
    //                     let base_qual = unsafe { *read.qual().get_unchecked(qpos as usize) };
    //                     // Prob_read_base: Wkeit, dass die gegebene Readbase tatsachlich der 2. base entspricht (Also dass es eigtl die 2. Base ist)
    //                     prob_alt = prob_read_base(read_base, b'C', base_qual);
    //                     let no_c = if read_base != b'C' { read_base } else { b'T' };
    //                     prob_ref = prob_read_base(read_base, no_c, base_qual);
    //                 }
    //                 else {
    //                     let read_base = unsafe { read.seq().decoded_base_unchecked((qpos + 1) as usize) };
    //                     let base_qual = unsafe { *read.qual().get_unchecked((qpos + 1) as usize) };
    //                     prob_alt = prob_read_base(read_base, b'G', base_qual);
    //                     let no_g = if read_base != b'G' { read_base } else { b'A' };
    //                     prob_ref = prob_read_base(read_base, no_g, base_qual);
    //                 }
    //             }
    //             Readtype::PacBio => {
    //                 // Get methylation info from MM and ML TAG.
    //                 let meth_pos = meth_pos(read).unwrap();
    //                 let meth_probs = meth_probs(read).unwrap();
    //                 let pos_to_probs: HashMap<usize, f64> =
    //                     meth_pos.into_iter().zip(meth_probs.into_iter()).collect();
    //                 if let Some(value) = pos_to_probs.get(&(qpos as usize)) {
    //                     prob_alt = LogProb::from(Prob(*value as f64));
    //                     prob_ref = LogProb::from(Prob(1 as f64 - *value as f64));
    //                 } else {
    //                     // TODO What should I do if there is no prob given
    //                     prob_alt = LogProb::from(Prob(0.0));
    //                     prob_ref = LogProb::from(Prob(1.0));
    //                     warn!("No probability given for unmethylated Cs!");
    //                 }
    //             }
    //         }
    //         // TODO: Implement strand
    //         let strand = if prob_ref != prob_alt {
    //             Strand::from_record_and_pos(read, qpos as usize)?
    //         } else {
    //             // METHOD: if record is not informative, we don't want to
    //             // retain its information (e.g. strand).
    //             Strand::no_strand_info()
    //         };
    //         // let strand = Strand::no_strand_info();
    //         Ok(Some(
    //             AlleleSupportBuilder::default()
    //                 .prob_ref_allele(prob_ref)
    //                 .prob_alt_allele(prob_alt)
    //                 .strand(strand)
    //                 .read_position(Some(qpos as u32))
    //                 // TODO: Implement third allele
    //                 .third_allele_evidence(None)
    //                 .build()
    //                 .unwrap(),
    //         ))
    //     } else {
    //         // a read that spans an SNV might have the respective position in the
    //         // reference skipped (Cigar op 'N'), and the library should not choke on those reads
    //         // but instead needs to know NOT to add those reads (as observations) further up
    //         Ok(None)
    //     }
    // }
    // #################################################################### - Old SingleEnd - ###################################################




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

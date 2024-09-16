// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use anyhow::Result;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};
use rust_htslib::bam;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::bases::prob_read_base_miscall;
use crate::variants::evidence::observations::read_observation::Strand;
use crate::variants::evidence::realignment::Realignable;
use crate::variants::model;
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, Evidence, Overlap, SingleLocus, Variant,
};

use super::{MultiLocus, ToVariantRepresentation};

#[derive(Debug)]
pub(crate) struct None {
    loci: MultiLocus,
    ref_base: u8,
}

impl None {
    pub(crate) fn new(locus: genome::Locus, ref_base: u8) -> Self {
        None {
            loci: MultiLocus::from_single_locus(SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + 1,
            ))),
            ref_base: ref_base.to_ascii_uppercase(),
        }
    }

    fn allele_support_per_read(&self, read: &bam::Record) -> Result<Option<AlleleSupport>> {
        if let Overlap::Enclosing = self.locus().overlap(read, false) {
            return Ok(None);
        }

        if let Some(qpos) = read
            .cigar_cached()
            .unwrap()
            // TODO expect u64 in read_pos
            .read_pos(self.locus().range().start as u32, false, false)?
        {
            let read_base = read.seq()[qpos as usize].to_ascii_uppercase();
            let base_qual = read.qual()[qpos as usize];
            let prob_miscall = prob_read_base_miscall(base_qual);

            let (prob_ref, prob_alt) = if read_base == self.ref_base {
                (prob_miscall.ln_one_minus_exp(), prob_miscall)
            } else {
                (prob_miscall, prob_miscall.ln_one_minus_exp())
            };

            let strand = if prob_ref != prob_alt {
                Strand::from_record_and_pos(read, qpos as usize)?
            } else {
                // METHOD: if record is not informative, we don't want to
                // retain its information (e.g. strand).
                Strand::None
            };

            Ok(Some(
                AlleleSupportBuilder::default()
                    .prob_ref_allele(prob_ref)
                    .prob_alt_allele(prob_alt)
                    .strand(strand)
                    .build()
                    .unwrap(),
            ))
        } else {
            // a read that spans a potential Ref site might have the respective position deleted (Cigar op 'D')
            // or reference skipped (Cigar op 'N'), and the library should not choke on those reads
            // but instead needs to know NOT to add those reads (as observations) further up
            Ok(None)
        }
    }

    fn locus(&self) -> &SingleLocus {
        &self.loci[0]
    }
}

impl Variant for None {
    fn is_imprecise(&self) -> bool {
        false
    }

    fn is_valid_evidence(
        &self,
        evidence: &Evidence,
        _: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        match evidence {
            Evidence::SingleEnd(read) => {
                if let Overlap::Enclosing = self.locus().overlap(read, false) {
                    Some(vec![0])
                } else {
                    None
                }
            }
            Evidence::PairedEnd { left, right } => {
                if let Overlap::Enclosing = self.locus().overlap(left, false) {
                    Some(vec![0])
                } else if let Overlap::Enclosing = self.locus().overlap(right, false) {
                    Some(vec![0])
                } else {
                    None
                }
            }
        }
    }

    fn loci(&self) -> &MultiLocus {
        &self.loci
    }

    fn allele_support(
        &self,
        evidence: &Evidence,
        _: &AlignmentProperties,
        _: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        match evidence {
            Evidence::SingleEnd(read) => Ok(self.allele_support_per_read(read)?),
            Evidence::PairedEnd { left, right } => {
                let left_support = self.allele_support_per_read(left)?;
                let right_support = self.allele_support_per_read(right)?;

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
        }
    }

    fn prob_sample_alt(&self, _: &Evidence, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

impl ToVariantRepresentation for None {
    fn to_variant_representation(&self) -> model::Variant {
        model::Variant::None
    }
}

// #[cfg(test)]
// mod tests {

//     use super::*;
//     use crate::model;

//     use rust_htslib::bam::record::{Cigar, CigarString};
//     use std::str;

//     #[test]
//     fn test_prob_none() {
//         let ref_seq: Vec<u8> = b"GATTACA"[..].to_owned();

//         let mut records: Vec<bam::Record> = Vec::new();
//         let mut qname: &[u8];
//         let mut seq: &[u8];

//         let mut none_evidence = NoneEvidence::new();

//         // Ignore leading HardClip, skip leading SoftClip, reference nucleotide
//         qname = b"HC_SC_ref";
//         let cigar = CigarString(vec![
//             Cigar::HardClip(3),
//             Cigar::SoftClip(1),
//             Cigar::Match(5),
//         ]);
//         seq = b"TATTaC";
//         let qual = [20, 30, 30, 30, 40, 30];
//         let mut record1 = bam::Record::new();
//         record1.set(qname, Some(&cigar), seq, &qual);
//         record1.set_pos(1);
//         records.push(record1);

//         // Ignore leading HardClip, skip leading SoftClip, non-reference nucleotide
//         qname = b"HC_SC_non-ref";
//         let cigar = CigarString(vec![
//             Cigar::HardClip(5),
//             Cigar::SoftClip(2),
//             Cigar::Match(4),
//         ]);
//         seq = b"TTTTCC";
//         let qual = [15, 15, 20, 20, 30, 20];
//         let mut record2 = bam::Record::new();
//         record2.set(qname, Some(&cigar), seq, &qual);
//         record2.set_pos(2);
//         records.push(record2);

//         // reference nucleotide, trailing SoftClip, trailing HardClip
//         qname = b"ref_SC_HC";
//         let cigar = CigarString(vec![
//             Cigar::Match(3),
//             Cigar::SoftClip(2),
//             Cigar::HardClip(7),
//         ]);
//         seq = b"ACATA";
//         let qual = [50, 20, 20, 20, 20];
//         let mut record3 = bam::Record::new();
//         record3.set(qname, Some(&cigar), seq, &qual);
//         record3.set_pos(4);
//         records.push(record3);

//         // three nucleotide Deletion covering Ref position
//         qname = b"M_3Del_M";
//         let cigar = CigarString(vec![Cigar::Match(3), Cigar::Del(3), Cigar::Match(1)]);
//         seq = b"GATA";
//         let qual = [10, 30, 30, 30];
//         let mut record4 = bam::Record::new();
//         record4.set(qname, Some(&cigar), seq, &qual);
//         record4.set_pos(0);
//         records.push(record4);

//         // truth
//         let probs_ref = [0.9999, 0.001, 0.99999];
//         let probs_alt = [0.0001, 0.999, 0.00001];
//         let eps = [0.00001, 0.0001, 0.000001];

//         let vpos = 4;
//         let variant = model::Variant::None;
//         for (i, mut rec) in records.into_iter().enumerate() {
//             rec.cache_cigar();
//             println!("{}", str::from_utf8(rec.qname()).unwrap());
//             if let Ok(Some((prob_ref, prob_alt))) =
//                 none_evidence.prob(&rec, rec.cigar_cached().unwrap(), vpos, &variant, &ref_seq)
//             {
//                 println!("{:?}", rec.cigar_cached());
//                 println!(
//                     "Pr(ref)={} Pr(alt)={}",
//                     (*prob_ref).exp(),
//                     (*prob_alt).exp()
//                 );
//                 assert_relative_eq!((*prob_ref).exp(), probs_ref[i], epsilon = eps[i]);
//                 assert_relative_eq!((*prob_alt).exp(), probs_alt[i], epsilon = eps[i]);
//             } else {
//                 // tests for reference position not being covered should be pushed onto records last
//                 // and should have 10 as the quality value of the first base in seq
//                 assert_eq!(rec.qual()[0], 10);
//             }
//         }
//     }
// }

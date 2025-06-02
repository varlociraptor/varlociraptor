// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::cmp;
use std::ops::Range;

use std::sync::Arc;

use anyhow::Result;

use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};
use rust_htslib::bam;

use crate::default_ref_base_emission;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils;
use crate::variants::evidence::bases::prob_read_base;
use crate::variants::evidence::observations::read_observation::Strand;
use crate::variants::evidence::realignment::edit_distance::EditDistance;
use crate::variants::evidence::realignment::pairhmm::RefBaseEmission;
use crate::variants::evidence::realignment::pairhmm::RefBaseVariantEmission;
use crate::variants::evidence::realignment::pairhmm::VariantEmission;
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::model;
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, Evidence, Overlap, SingleLocus, Variant,
};

use super::MultiLocus;
use super::ToVariantRepresentation;

#[derive(Debug)]
pub(crate) struct Snv<R: Realigner> {
    loci: MultiLocus,
    ref_base: u8,
    alt_base: u8,
    realigner: RefCell<R>,
    realign_indel_reads: bool,
}

impl<R: Realigner> Snv<R> {
    pub(crate) fn new(
        locus: genome::Locus,
        ref_base: u8,
        alt_base: u8,
        realigner: R,
        realign_indel_reads: bool,
    ) -> Self {
        Snv {
            loci: MultiLocus::from_single_locus(SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + 1,
            ))),
            ref_base: ref_base.to_ascii_uppercase(),
            alt_base: alt_base.to_ascii_uppercase(),
            realigner: RefCell::new(realigner),
            realign_indel_reads,
        }
    }

    fn allele_support_per_read(
        &self,
        read: &bam::Record,
        alignment_properties: &AlignmentProperties,
        alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        if self.locus().overlap(read, false, 0, 0) != Overlap::Enclosing {
            return Ok(None);
        }

        if self.realign_indel_reads && utils::contains_indel_op(read) {
            // METHOD: reads containing indel operations should always be realigned,
            // as their support or non-support of the SNV might be an artifact
            // of the aligner.
            Ok(Some(self.realigner.borrow_mut().allele_support(
                read,
                self.loci.iter(),
                self,
                alt_variants,
                alignment_properties,
            )?))
        } else if let Some(qpos) = read
            .cigar_cached()
            .unwrap()
            // TODO expect u64 in read_pos
            .read_pos(self.locus().range().start as u32, false, false)?
        {
            let read_base =
                unsafe { read.seq().decoded_base_unchecked(qpos as usize) }.to_ascii_uppercase();
            let base_qual = unsafe { *read.qual().get_unchecked(qpos as usize) };
            let prob_alt = prob_read_base(read_base, self.alt_base, base_qual);
            let mut is_third_allele = false;

            // METHOD: instead of considering the actual REF base, we assume that REF is whatever
            // base the read has at this position (if not the ALT base). This way, we avoid biased
            // allele frequencies at sites with multiple alternative alleles.
            // Note that this is an approximation. The real solution would be to have multiple allele
            // frequency variables in the likelihood function, but that would be computationally
            // more demanding (leading to a combinatorial explosion).
            // However, the approximation is pretty accurate, because it will only matter for true
            // multiallelic cases. Sequencing errors won't have a severe effect on the allele frequencies
            // because they are too rare.
            // Here, N bases do not count as additional edits that would indicate a third allele.
            let non_alt_base = if read_base != b'N' && read_base != self.alt_base {
                is_third_allele = read_base != self.ref_base;
                read_base
            } else {
                self.ref_base
            };

            let prob_ref = prob_read_base(read_base, non_alt_base, base_qual);
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
                    // METHOD: hardclips are not part of qpos, but they are part
                    // of the original read sequence. Hence, they have to be added
                    // here.
                    .read_position(Some(
                        qpos + read.cigar_cached().unwrap().leading_hardclips() as u32,
                    ))
                    .third_allele_evidence(if is_third_allele {
                        Some(EditDistance(1))
                    } else {
                        None
                    })
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

    fn locus(&self) -> &SingleLocus {
        &self.loci[0]
    }
}

impl<R: Realigner> Realignable for Snv<R> {
    fn alt_emission_params(
        &self,
        ref_buffer: Arc<reference::Buffer>,
        _: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<Box<dyn RefBaseVariantEmission>>> {
        let start = self.locus().range().start as usize;

        let ref_seq = ref_buffer.seq(self.locus().contig())?;

        let ref_seq_len = ref_seq.len();
        Ok(vec![Box::new(SnvEmissionParams {
            ref_seq,
            ref_offset: start.saturating_sub(ref_window),
            ref_end: cmp::min(start + 1 + ref_window, ref_seq_len),
            alt_start: start,
            alt_base: self.alt_base,
            ref_offset_override: None,
            ref_end_override: None,
        })])
    }
}

impl<R: Realigner> Variant for Snv<R> {
    fn is_imprecise(&self) -> bool {
        false
    }

    fn is_valid_evidence(
        &self,
        evidence: &Evidence,
        _: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        match evidence {
            Evidence::SingleEndSequencingRead(read) => {
                if let Overlap::Enclosing = self.locus().overlap(read, false, 0, 0) {
                    Some(vec![0])
                } else {
                    None
                }
            }
            Evidence::PairedEndSequencingRead { left, right } => {
                if let Overlap::Enclosing = self.locus().overlap(left, false, 0, 0) {
                    Some(vec![0])
                } else if let Overlap::Enclosing = self.locus().overlap(right, false, 0, 0) {
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
        alignment_properties: &AlignmentProperties,
        alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        match evidence {
            Evidence::SingleEndSequencingRead(read) => {
                Ok(self.allele_support_per_read(read, alignment_properties, alt_variants)?)
            }
            Evidence::PairedEndSequencingRead { left, right } => {
                let left_support =
                    self.allele_support_per_read(left, alignment_properties, alt_variants)?;
                let right_support =
                    self.allele_support_per_read(right, alignment_properties, alt_variants)?;

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

impl<R: Realigner> ToVariantRepresentation for Snv<R> {
    fn to_variant_representation(&self) -> model::Variant {
        model::Variant::Snv(self.alt_base)
    }
}

/// Emission parameters for PairHMM over insertion allele.
pub(crate) struct SnvEmissionParams {
    ref_seq: Arc<Vec<u8>>,
    ref_offset: usize,
    ref_end: usize,
    alt_start: usize,
    alt_base: u8,
    ref_offset_override: Option<usize>,
    ref_end_override: Option<usize>,
}

impl RefBaseEmission for SnvEmissionParams {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;

        if i_ != self.alt_start {
            self.ref_seq[i_]
        } else {
            self.alt_base
        }
    }

    fn variant_homopolymer_ref_range(&self) -> Option<Range<u64>> {
        None
    }

    fn variant_ref_range(&self) -> Option<Range<u64>> {
        Some(self.alt_start as u64..self.alt_start as u64 + 1)
    }

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }

    default_ref_base_emission!();
}

impl VariantEmission for SnvEmissionParams {
    fn is_homopolymer_indel(&self) -> bool {
        false
    }

    fn alt_vs_ref_len_diff(&self) -> isize {
        0
    }
}

// #[cfg(test)]
// mod tests {

//     use super::*;
//     use crate::model;

//     use rust_htslib::bam;
//     use rust_htslib::bam::record::{Cigar, CigarString};
//     use std::str;

//     #[test]
//     fn test_prob_snv() {
//         let ref_seq: Vec<u8> = b"CCTATACGCGT"[..].to_owned();

//         let mut records: Vec<bam::Record> = Vec::new();
//         let mut qname: &[u8];
//         let mut seq: &[u8];

//         let mut snv_evidence = SNVEvidence::new();

//         // Ignore leading HardClip, skip leading SoftClip, reference nucleotide
//         qname = b"HC_SC_M";
//         let cigar = CigarString(vec![
//             Cigar::HardClip(5),
//             Cigar::SoftClip(2),
//             Cigar::Match(6),
//         ]);
//         seq = b"AATATACG";
//         let qual = [20, 20, 30, 30, 30, 40, 30, 30];
//         let mut record1 = bam::Record::new();
//         record1.set(qname, Some(&cigar), seq, &qual);
//         record1.set_pos(2);
//         records.push(record1);

//         // Ignore leading HardClip, skip leading Insertion, alternative nucleotide
//         qname = b"HC_Ins_M";
//         let cigar = CigarString(vec![Cigar::HardClip(2), Cigar::Ins(2), Cigar::Match(6)]);
//         seq = b"TTTATGCG";
//         let qual = [20, 20, 20, 20, 20, 30, 20, 20];
//         let mut record2 = bam::Record::new();
//         record2.set(qname, Some(&cigar), seq, &qual);
//         record2.set_pos(2);
//         records.push(record2);

//         // Matches and deletion before position, reference nucleotide
//         qname = b"Eq_Diff_Del_Eq";
//         let cigar = CigarString(vec![
//             Cigar::Equal(2),
//             Cigar::Diff(1),
//             Cigar::Del(2),
//             Cigar::Equal(5),
//         ]);
//         seq = b"CCAACGCG";
//         let qual = [30, 30, 30, 50, 30, 30, 30, 30];
//         let mut record3 = bam::Record::new();
//         record3.set(qname, Some(&cigar), seq, &qual);
//         record3.set_pos(0);
//         records.push(record3);

//         // single nucleotide Deletion covering SNV position
//         qname = b"M_Del_M";
//         let cigar = CigarString(vec![Cigar::Match(4), Cigar::Del(1), Cigar::Match(4)]);
//         seq = b"CTATCGCG";
//         let qual = [10, 30, 30, 30, 30, 30, 30, 30];
//         let mut record4 = bam::Record::new();
//         record4.set(qname, Some(&cigar), seq, &qual);
//         record4.set_pos(1);
//         records.push(record4);

//         // three nucleotide RefSkip covering SNV position
//         qname = b"M_RefSkip_M";
//         let cigar = CigarString(vec![
//             Cigar::Equal(1),
//             Cigar::Diff(1),
//             Cigar::Equal(2),
//             Cigar::RefSkip(3),
//             Cigar::Match(4),
//         ]);
//         seq = b"CTTAGCGT";
//         let qual = [10, 30, 30, 30, 30, 30, 30, 30];
//         let mut record5 = bam::Record::new();
//         record5.set(qname, Some(&cigar), seq, &qual);
//         record5.set_pos(0);
//         records.push(record5);

//         // truth
//         let probs_ref = [0.9999, 0.00033, 0.99999];
//         let probs_alt = [0.000033, 0.999, 0.0000033];
//         let eps = [0.000001, 0.00001, 0.0000001];

//         let vpos = 5;
//         let variant = model::Variant::SNV(b'G');
//         for (i, mut rec) in records.into_iter().enumerate() {
//             rec.cache_cigar();
//             println!("{}", str::from_utf8(rec.qname()).unwrap());
//             if let Ok(Some((prob_ref, prob_alt))) =
//                 snv_evidence.prob(&rec, rec.cigar_cached().unwrap(), vpos, &variant, &ref_seq)
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

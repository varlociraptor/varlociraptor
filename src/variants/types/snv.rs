// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::cmp;
use std::rc::Rc;
use std::sync::Arc;

use anyhow::Result;
use bio::stats::pairhmm::EmissionParameters;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils;
use crate::variants::evidence::bases::prob_read_base;
use crate::variants::evidence::observation::Strand;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, Overlap, SingleEndEvidence, SingleLocus, Variant,
};
use crate::{default_emission, default_ref_base_emission};

pub(crate) struct Snv<R: Realigner> {
    locus: SingleLocus,
    ref_base: u8,
    alt_base: u8,
    realigner: RefCell<R>,
}

impl<R: Realigner> Snv<R> {
    pub(crate) fn new(locus: genome::Locus, ref_base: u8, alt_base: u8, realigner: R) -> Self {
        Snv {
            locus: SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + 1,
            )),
            ref_base: ref_base.to_ascii_uppercase(),
            alt_base: alt_base.to_ascii_uppercase(),
            realigner: RefCell::new(realigner),
        }
    }
}

impl<'a, R: Realigner> Realignable<'a> for Snv<R> {
    type EmissionParams = SnvEmissionParams<'a>;

    fn alt_emission_params(
        &self,
        read_emission_params: Rc<ReadEmission<'a>>,
        ref_buffer: Arc<reference::Buffer>,
        _: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<SnvEmissionParams<'a>>> {
        let start = self.locus.range().start as usize;

        let ref_seq = ref_buffer.seq(self.locus.contig())?;

        let ref_seq_len = ref_seq.len();
        Ok(vec![SnvEmissionParams {
            ref_seq,
            ref_offset: start.saturating_sub(ref_window),
            ref_end: cmp::min(start + 1 + ref_window, ref_seq_len),
            alt_start: start,
            alt_base: self.alt_base,
            read_emission: read_emission_params,
        }])
    }
}

impl<R: Realigner> Variant for Snv<R> {
    type Evidence = SingleEndEvidence;
    type Loci = SingleLocus;

    fn is_valid_evidence(
        &self,
        evidence: &SingleEndEvidence,
        _: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        if let Overlap::Enclosing = self.locus.overlap(evidence, false) {
            Some(vec![0])
        } else {
            None
        }
    }

    fn loci(&self) -> &SingleLocus {
        &self.locus
    }

    fn allele_support(
        &self,
        read: &SingleEndEvidence,
        _: &AlignmentProperties,
    ) -> Result<Option<AlleleSupport>> {
        if utils::contains_indel_op(&**read) {
            // METHOD: reads containing indel operations should always be realigned,
            // as their support or non-support of the SNV might be an artifact
            // of the aligner.
            Ok(Some(self.realigner.borrow_mut().allele_support(
                &**read,
                [&self.locus].iter(),
                self,
            )?))
        } else if let Some(qpos) = read
            .cigar_cached()
            .unwrap()
            // TODO expect u64 in read_pos
            .read_pos(self.locus.range().start as u32, false, false)?
        {
            let read_base = unsafe { read.seq().decoded_base_unchecked(qpos as usize) };
            let base_qual = unsafe { *read.qual().get_unchecked(qpos as usize) };
            let prob_alt = prob_read_base(read_base, self.alt_base, base_qual);

            // METHOD: instead of considering the actual REF base, we assume that REF is whatever
            // base the read has at this position (if not the ALT base). This way, we avoid biased
            // allele frequencies at sites with multiple alternative alleles.
            // Note that this is an approximation. The real solution would be to have multiple allele
            // frequency variables in the likelihood function, but that would be computationally
            // more demanding (leading to a combinatorial explosion).
            // However, the approximation is pretty accurate, because it will only matter for true
            // multiallelic cases. Sequencing errors won't have a severe effect on the allele frequencies
            // because they are too rare.
            let non_alt_base = if read_base != self.alt_base {
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
                    .read_position(Some(qpos))
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

    fn prob_sample_alt(&self, _: &SingleEndEvidence, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

/// Emission parameters for PairHMM over insertion allele.
pub(crate) struct SnvEmissionParams<'a> {
    ref_seq: Arc<Vec<u8>>,
    ref_offset: usize,
    ref_end: usize,
    alt_start: usize,
    alt_base: u8,
    read_emission: Rc<ReadEmission<'a>>,
}

impl<'a> RefBaseEmission for SnvEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;

        if i_ != self.alt_start {
            self.ref_seq[i_]
        } else {
            self.alt_base
        }
    }

    default_ref_base_emission!();
}

impl<'a> EmissionParameters for SnvEmissionParams<'a> {
    default_emission!();

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
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

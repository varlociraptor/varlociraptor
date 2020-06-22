use std::str;

use anyhow::Result;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::utils::is_reverse_strand;
use crate::variants::evidence::bases::prob_read_base;
use crate::variants::{
    AlleleSupport, AlleleSupportBuilder, Overlap, SingleEndEvidence, SingleLocus, Variant,
};

pub struct SNV {
    locus: SingleLocus,
    ref_base: u8,
    alt_base: u8,
}

impl SNV {
    pub fn new(locus: genome::Locus, ref_base: u8, alt_base: u8) -> Self {
        SNV {
            locus: SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + 1,
            )),
            ref_base: ref_base.to_ascii_uppercase(),
            alt_base: alt_base.to_ascii_uppercase(),
        }
    }
}

impl Variant for SNV {
    type Evidence = SingleEndEvidence;
    type Loci = SingleLocus;

    fn is_valid_evidence(&self, evidence: &SingleEndEvidence) -> Option<Vec<usize>> {
        dbg!((str::from_utf8(evidence.qname()).unwrap(), self.locus.overlap(evidence, false)));
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
        if let Some(qpos) = read
            .cigar_cached()
            .unwrap()
            // TODO expect u64 in read_pos
            .read_pos(self.locus.range().start as u32, false, false)?
        {
            let read_base = read.seq()[qpos as usize].to_ascii_uppercase();
            let base_qual = read.qual()[qpos as usize];
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

            Ok(Some(
                AlleleSupportBuilder::default()
                    .prob_ref_allele(prob_ref)
                    .prob_alt_allele(prob_alt)
                    .register_record(read)
                    .build()
                    .unwrap(),
            ))
        } else {
            // a read that spans an SNV might have the respective position deleted (Cigar op 'D')
            // or reference skipped (Cigar op 'N'), and the library should not choke on those reads
            // but instead needs to know NOT to add those reads (as observations) further up
            Ok(None)
        }
    }

    fn prob_sample_alt(&self, _: &SingleEndEvidence, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
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

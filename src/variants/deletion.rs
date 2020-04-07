use anyhow::Result;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::{AlleleProb, PairedEndEvidence, MultiLocus, SingleLocus, Variant};

pub struct Deletion {
    locus: genome::Interval,
    fetch_loci: MultiLocus,
}

impl Deletion {
    pub fn new(locus: genome::Interval) -> Self {
        let start = locus.range().start;
        let end = locus.range().end;
        let centerpoint = start + ((end - start) as f64 / 2.0).round() as u64;
        let contig = locus.contig().to_owned();

        let fetch_loci = MultiLocus::new(vec![
            SingleLocus(genome::Interval::new(contig.clone(), start..start + 1)), 
            SingleLocus(genome::Interval::new(contig.clone(), centerpoint, centerpoint + 1), 
            SingleLocus(genome::Interval::new(contig, end - 1..end)))
        ]);

        Deletion {
            locus: locus,
            fetch_loci,
        }
    }
}


impl<'a> Variant<'a> for Deletion {
    type Evidence = PairedEndEvidence<'a>;
    type Loci = MultiLocus;

    fn is_valid_evidence(&self, evidence: &Self::Evidence) -> Option<Vec<usize>> {
        let mut locus_idx = Vec::new();
        for locus in self.loci() {
            match evidence {
                PairedEndEvidence::SingleEnd(read) => {

                },
                PairedEndEvidence::PairedEnd { left, right } => {
                    
                }
            }
            let enclosing = pos < start && end_pos > end;
            if enclosing {
                Overlap::Enclosing(l)
            } else {
                if end_pos <= end && end_pos > start {
                    Overlap::Right(end_pos - start)
                } else if pos >= start && pos < end {
                    Overlap::Left(end - pos)
                } else {
                    Overlap::None
                }
            }
        }
        if locus_idx.is_empty() {
            None
        } else {
            Some(locus_idx)
        }
    }

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci {
        &self.fetch_loci
    }

    /// Calculate probability for alt and reference allele.
    fn prob_alleles(&self, evidence: &Self::Evidence) -> Result<Option<AlleleProb>> {

    }
}
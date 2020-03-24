use anyhow::Result;
use bio::stats::{LogProb, PHREDProb, Prob};
use rust_htslib::bam;

use crate::utils::{Overlap, GenomicLocus};
use crate::model::sample;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::model::evidence::observation::{Observable, ObservationBuilder, Observation, Evidence, SingleEndEvidence, PairedEndEvidence, Strand, StrandBuilder};

pub mod mnv;
pub mod snv;

#[derive(Debug, CopyGetters, new)]
#[getset(get_copy = "pub")]
pub struct AlleleProb {
    ref_allele: LogProb,
    alt_allele: LogProb,
}

impl AlleleProb {
    pub fn missed_allele(&self) -> LogProb {
        self.ref_allele.ln_add_exp(self.alt_allele) - LogProb(2.0_f64.ln())
    }
}

pub trait Variant<'a> {
    type Evidence: Evidence<'a>;
    type Loci: Loci;

    /// Determine whether record overlaps the variant.
    fn overlap(&self, evidence: &Self::Evidence) -> Overlap;

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci;

    /// Calculate probability for alt and reference allele.
    fn prob_alleles(&self, evidence: &Self::Evidence) -> Result<Option<AlleleProb>>;

    /// Calculate probability to sample a record length like the given one from the alt allele.
    fn prob_sample_alt(&self, evidence: &Self::Evidence, alignment_properties: &AlignmentProperties) -> LogProb;
}

impl<'a, V> Observable<'a, SingleEndEvidence<'a>> for V where V: Variant<'a, Evidence=SingleEndEvidence<'a>, Loci=SingleLocus> {

    fn prob_mapping(&self, evidence: &SingleEndEvidence<'a>) -> LogProb {
        let prob_mismapping = LogProb::from(PHREDProb(evidence.mapq() as f64));
        let prob_mapping = prob_mismapping.ln_one_minus_exp();
        prob_mapping
    }
    
    fn strand(&self, evidence: &SingleEndEvidence<'a>) -> Strand {
        let reverse = evidence.flags() & 0x10 != 0;

        StrandBuilder::default()
            .reverse(reverse)
            .forward(!reverse)
            .build()
            .unwrap()
    }

    fn extract_observations(&self, buffer: &'a mut sample::RecordBuffer, alignment_properties: &AlignmentProperties, max_depth: usize) -> Result<Vec<Observation>> {
        let locus = self.loci();
        buffer.fetch_surrounding(locus.chrom(), locus.pos(), locus.pos() + locus.len())?;

        let candidates: Vec<_> = buffer.iter().filter_map(|record| {
            let evidence = SingleEndEvidence::new(record);
            if record.pos() as u32 > locus.pos() || !self.overlap(&evidence).is_enclosing() {
                None
            } else {
                Some(evidence)
            }
        }).collect();

        let mut subsample = sample::SubsampleCandidates::new(max_depth, candidates.len());

        candidates.iter().filter_map(|evidence| {
            if subsample.keep() {
                match self.prob_alleles(evidence) {
                    Ok(Some(prob_alleles)) => {
                        let strand = self.strand(evidence);

                        Some(Ok(ObservationBuilder::default()
                            .prob_mapping_mismapping(self.prob_mapping(evidence))
                            .prob_alt(prob_alleles.alt_allele())
                            .prob_ref(prob_alleles.ref_allele())
                            .prob_sample_alt(self.prob_sample_alt(evidence, alignment_properties))
                            .prob_missed_allele(prob_alleles.missed_allele())
                            .prob_overlap(LogProb::ln_zero()) // no double overlap possible
                            .prob_any_strand(LogProb::from(Prob(0.5)))
                            .forward_strand(strand.forward())
                            .reverse_strand(strand.reverse())
                            .build()
                            .unwrap()))
                    },
                    Ok(None) => None,
                    Err(e) => Some(Err(e))
                }
            } else {
                None
            }
        }).collect::<Result<Vec<Observation>>>()
    }
}

impl<'a, V> Observable<'a, PairedEndEvidence<'a>> for V where V: Variant<'a, Evidence=PairedEndEvidence<'a>, Loci=SingleLocus> {

    fn prob_mapping(&self, evidence: &PairedEndEvidence<'a>) -> LogProb {
        let prob = |record: &bam::Record| LogProb::from(PHREDProb(record.mapq() as f64));
        match evidence {
            PairedEndEvidence::SingleEnd(record) => prob(record),
            PairedEndEvidence::PairedEnd{ left, right } => prob(left) + prob(right),
        }
    }
    
    fn strand(&self, evidence: &PairedEndEvidence<'a>) -> Strand {
        let reverse = evidence.flags() & 0x10 != 0;

        StrandBuilder::default()
            .reverse(reverse)
            .forward(!reverse)
            .build()
            .unwrap()
    }
}


pub trait Loci {}

#[derive(new, Debug, Getters, CopyGetters, Derefable)]
pub struct SingleLocus {
    #[getset(get = "pub")]
    #[deref]
    locus: GenomicLocus,
    #[getset(get_copy = "pub")]
    len: u32
}

impl Loci for SingleLocus {}

#[derive(new, Debug)]
pub struct MultiLocus {
    loci: Vec<SingleLocus>,
}

impl Loci for MultiLocus {}
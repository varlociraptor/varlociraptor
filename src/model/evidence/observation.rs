use std::str;
use std::rc::Rc;
use std::ops::Range;

use vec_map::VecMap;
use itertools::Itertools;
use rgsl::randist::poisson::poisson_pdf;
use serde::Serialize;
use serde::ser::{Serializer, SerializeStruct};

use bio::stats::LogProb;
use rust_htslib::bam;
use rust_htslib::bam::record::CigarString;

use model::{AlleleFreq, Variant};
use model::sample::Sample;


/// An observation for or against a variant.
#[derive(Clone)]
pub struct Observation<'a> {
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    pub prob_mapping: LogProb,
    /// Probability that the read/read-pair comes from the alternative allele.
    pub prob_alt: LogProb,
    /// Probability that the read/read-pair comes from the reference allele.
    pub prob_ref: LogProb,
    /// Probability of the read/read-pair given that it has been mismapped.
    pub prob_mismapped: LogProb,
    /// Read lengths
    pub left_read_len: u32,
    pub right_read_len: Option<u32>,
    /// Common stuff shared between observations
    pub common: Rc<Common<'a>>,
    /// Type of evidence.
    pub evidence: Evidence
}


impl<'a> Observation<'a> {
    pub fn is_alignment_evidence(&self) -> bool {
        if let Evidence::Alignment(_) = self.evidence {
            true
        } else {
            false
        }
    }

    pub fn prob_sample_alt(&self, allele_freq: AlleleFreq) -> LogProb {
        self.common.prob_sample_alt(
            allele_freq,
            self.left_read_len,
            self.right_read_len
        )
    }
}


impl<'a> Serialize for Observation<'a> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: Serializer {
        let mut s = serializer.serialize_struct("Observation", 3)?;
        s.serialize_field("prob_mapping", &self.prob_mapping)?;
        s.serialize_field("prob_alt", &self.prob_alt)?;
        s.serialize_field("prob_ref", &self.prob_ref)?;
        s.serialize_field("evidence", &self.evidence)?;
        s.end()
    }
}


/// Types of evidence that lead to an observation.
/// The contained information is intended for debugging and will be printed together with
/// observations.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Evidence {
    /// Insert size of fragment
    InsertSize(String),
    /// Alignment of a single read
    Alignment(String)
}


impl Evidence {
    /// Create a dummy alignment.
    pub fn dummy_alignment() -> Self {
        Evidence::Alignment("Dummy-Alignment".to_owned())
    }

    /// Create dummy insert size evidence.
    pub fn dummy_insert_size(insert_size: u32) -> Self {
        Evidence::InsertSize(format!("insert-size={}", insert_size))
    }

    /// Create insert size evidence.
    pub fn insert_size(
        insert_size: u32,
        left: &CigarString,
        right: &CigarString,
        left_record: &bam::Record,
        right_record: &bam::Record,
        p_left_ref: LogProb,
        p_left_alt: LogProb,
        p_right_ref: LogProb,
        p_right_alt: LogProb,
        p_isize_ref: LogProb,
        p_isize_alt: LogProb
    ) -> Self {
        Evidence::InsertSize(format!(
            "left: cigar={} ({:e} vs {:e}), right: cigar={} ({:e} vs {:e}), insert-size={} ({:e} vs {:e}), qname={}, left: AS={:?}, XS={:?}, right: AS={:?}, XS={:?}",
            left, p_left_ref.exp(), p_left_alt.exp(),
            right, p_right_ref.exp(), p_right_alt.exp(),
            insert_size, p_isize_ref.exp(), p_isize_alt.exp(),
            str::from_utf8(left_record.qname()).unwrap(),
            left_record.aux(b"AS").map(|a| a.integer()),
            left_record.aux(b"XS").map(|a| a.integer()),
            right_record.aux(b"AS").map(|a| a.integer()),
            right_record.aux(b"XS").map(|a| a.integer())
        ))
    }

    /// Create alignment evidence.
    pub fn alignment(cigar: &CigarString, record: &bam::Record) -> Self {
        Evidence::Alignment(format!(
            "cigar={}, qname={}, AS={:?}, XS={:?}",
            cigar, str::from_utf8(record.qname()).unwrap(),
            record.aux(b"AS").map(|a| a.integer()),
            record.aux(b"XS").map(|a| a.integer())
        ))
    }
}


/// Data that is shared among all observations over a locus.
#[derive(Clone)]
pub struct Common<'a> {
    pub softclip_obs: VecMap<u32>,
    /// Average number of reads starting at any position in the region.
    pub coverage: f64,
    pub max_read_len: u32,
    pub enclosing_possible: bool,
    pub sample: &'a Sample,
    pub variant: Rc<Variant>
}


impl<'a> Common<'a> {
    pub fn new(
        sample: &'a Sample,
        variant: Rc<Variant>
    ) -> Self {
        Common {
            sample: sample,
            variant: variant,
            enclosing_possible: true,
            max_read_len: 0,
            coverage: 0.0,
            softclip_obs: VecMap::new()
        }
    }

    /// Calculate probability to sample reads from alt allele.
    /// For SNVs this is always 1.0.
    /// Otherwise, this has two components.
    ///
    /// # Component 1: maximum softclip
    /// First, we estimate the probability distribution of
    /// the maximum possible softclip. This is influenced by the implementation of the mapper
    /// and the repeat structure in the region of the variant. E.g., if a region is repetetive,
    /// a fragment from the alt allele that could map over the variant given the mapper sensitivity
    /// could still be placed somewhere else because there is a repeat that looks like the alt
    /// allele. Hence, the maximum softclip has to be calculated per region.
    ///
    /// Simply taking the maximum observed softclip in a region is not robust enough though,
    /// because small allele frequencies and low coverage can lead to not enough observations.
    /// Instead, we calculate the probability distribution of the maximum softclip, by assuming
    /// that read start positions are poisson distributed with mean given by the expected variant
    /// allele coverage.
    /// This coverage is calculated by multiplying allele frequency with the observed average
    /// number of reads starting at any position in the region.
    /// Then, we can calculate the likelihood of the observed softclips given a true maximum
    /// softclip by taking the product of the poisson distributed probabilities for each observed
    /// softclip count.
    /// By applying Bayes theorem, we obtain a posterior probability for each possible maximum
    /// softclip.
    ///
    /// # Component 2:
    /// We calculate the probability to sample a fragment from the alt allele given a maximum
    /// softclip and the read lengths. If the variant is small enough to be encoded in the CIGAR
    /// string, we can simply ignore the maximum softclip distribution.
    ///
    /// # Total probability
    /// The final result is obtained by combining the two components to a total probability.
    pub fn prob_sample_alt(
        &self,
        allele_freq: AlleleFreq,
        left_read_len: u32,
        right_read_len: Option<u32>
    ) -> LogProb {
        if self.variant.is_snv() {
            // For SNVs, sampling can be assumed to be always unbiased.
            return LogProb::ln_one();
        }

        let prob_sample_alt = |max_softclip| {
            if let Some(right_read_len) = right_read_len {
                // we have a proper pair, use fragement evidence
                self.sample.indel_fragment_evidence.borrow().prob_sample_alt(
                    left_read_len,
                    right_read_len,
                    max_softclip,
                    self.variant.as_ref()
                )
            } else {
                self.sample.indel_read_evidence.borrow().prob_sample_alt(
                    left_read_len,
                    max_softclip,
                    self.variant.as_ref()
                )
            }
        };

        if self.enclosing_possible {
            // if read can enclose variant, the maximum overlap is given by max read len
            prob_sample_alt(self.max_read_len)
        } else {
            // calculate total probability to sample alt allele given the max_softclip distribution

            // max_softclip likelihoods
            let likelihoods = self.softclip_range().map(
                |s| self.likelihood_max_softclip(s, allele_freq)
            ).collect_vec();
            let marginal = LogProb::ln_sum_exp(&likelihoods);

            LogProb::ln_sum_exp(&self.softclip_range().map(|max_softclip| {
                // posterior probability for maximum softclip s
                let prob_max_softclip = likelihoods[max_softclip as usize] - marginal;
                // probability to sample from alt allele with this maximum softclip
                prob_max_softclip + prob_sample_alt(max_softclip)
            }).collect_vec())
        }
    }

    fn softclip_range(&self) -> Range<u32> {
        0..self.max_read_len + 1
    }

    fn likelihood_max_softclip(
        &self,
        max_softclip: u32,
        allele_freq: AlleleFreq
    ) -> LogProb {
        let varcov = self.coverage * *allele_freq;
        self.softclip_range().map(|s| {
            let count = self.softclip_obs.get(s as usize).cloned().unwrap_or(0);
            let mu = if s <= max_softclip { varcov } else { 0.0 };
            LogProb(poisson_pdf(count, mu).ln())
        }).sum()
    }
}


impl<'a> Serialize for Common<'a> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: Serializer {
        let mut s = serializer.serialize_struct("CommonObservation", 2)?;
        s.serialize_field("softclip_obs", &self.softclip_obs)?;
        s.serialize_field("coverage", &self.coverage)?;
        s.end()
    }
}

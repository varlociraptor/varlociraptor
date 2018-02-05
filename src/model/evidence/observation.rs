use std::str;

use vec_map::VecMap;
use rgsl::randist::poisson::poisson_pdf;

use bio::stats::LogProb;
use rust_htslib::bam;
use rust_htslib::bam::record::CigarString;


/// An observation for or against a variant.
#[derive(Clone, Debug, Serialize, Deserialize)]
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
    pub right_read_len: u32,
    /// Common stuff shared between observations
    pub common: &'a Common,
    /// Type of evidence.
    pub evidence: Evidence
}


impl Observation {
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
pub struct Common<'a> {
    softclip_obs: VecMap<u32>,
    /// Average number of reads starting at any position in the region.
    coverage: f64,
    max_read_len: u32,
    enclosing_possible: bool,
    sample: &'a Sample,
    variant: &'a Variant
}


impl Common {
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
        right_read_len: u32
    ) -> LogProb {
        if variant.is_snv() {
            // For SNVs, sampling can be assumed to be always unbiased.
            return LogProb::ln_one();
        }

        if self.enclosing_possible {
            self.sample.indel_fragment_evidence.borrow().prob_sample_alt(
                left_read_len,
                right_read_len,
                // if read can enclose variant, the maximum overlap is given by max read len
                self.max_read_len,
                self.variant
            )
        } else {
            // calculate total probability to sample alt allele given the max_softclip distribution

            // max_softclip likelihoods
            let likelihoods = self.softclip_range().map(
                |s| self.likelihood_max_softclip(s, allele_freq)
            ).collect_vec();
            let marginal = LogProb::ln_sum_exp(likelihoods);

            LogProb::ln_sum_exp(&self.softclip_range().map(|max_softclip| {
                // posterior probability for maximum softclip s
                let prob_max_softclip = likelihoods[s as usize] - marginal;
                // probability to sample from alt allele given that softclip
                let prob_sample_alt = self.sample.indel_fragment_evidence.borrow().prob_sample_alt(
                    left_read_len,
                    right_read_len,
                    max_softclip,
                    self.variant
                );
                // joint probability
                prob_max_softclip + prob_sample_alt
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
        let varcov = self.coverage * allele_freq;
        let p = self.softclip_range().map(|s| {
            let count = self.softclip_obs.get(s).unwrap_or(0);
            let mu = if s <= max_softclip { varcov } else { 0 };
            LogProb(poisson_pdf(count, mu).ln())
        }).sum();
    }
}

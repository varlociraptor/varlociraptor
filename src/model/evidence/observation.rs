use std::str;
use std::collections::vec_deque;
use std::error::Error;
use std::cmp;
use std::rc::Rc;
use std::cell::RefCell;
use std::collections::BTreeMap;

use vec_map::VecMap;
use serde::Serialize;
use serde::ser::{Serializer, SerializeStruct};
use itertools::Itertools;
use rgsl::randist::poisson::poisson_pdf;

use bio::stats::LogProb;
use rust_htslib::bam;
use rust_htslib::bam::record::CigarString;

use model::Variant;
use model::evidence;
use model::sample::RecordBuffer;
use model::AlleleFreq;
use utils::{NUMERICAL_EPSILON, Overlap};


#[derive(Clone, Debug)]
pub enum ProbSampleAlt {
    /// probability depends on maximum possible softlip
    /// The given list of probabilities are the sampling probabilities given maximum
    /// softclips 0..probs.len().
    Dependent(Vec<LogProb>),
    /// probability does not depend on softclip
    Independent(LogProb),
    /// probability is always 1
    One
}


impl ProbSampleAlt {

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
    /// that read start positions are poisson distributed with a certain mean.
    /// The mean is calculated by counting start positions of reads with leading or trailing
    /// softclips that overlap the variant and dividing by the interval length between the
    /// first and the last start position. This is an estimate for the average number of softclipped
    /// reads from the variant allele per position.
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
    ///
    /// # Arguments
    /// * allele_freq - Given allele frequency of sample.
    pub fn joint_prob(&self, common_obs: &Common) -> LogProb {
        match self {
            &ProbSampleAlt::One => LogProb::ln_one(),
            &ProbSampleAlt::Independent(p) => p,
            &ProbSampleAlt::Dependent(ref probs) => {
                if common_obs.not_enough_softclips() {
                    // if there are no softclip observations we cannot infer anything
                    return LogProb::ln_one();
                }

                let softclip_range = || 0..probs.len() as u32;
                let varcov = common_obs.softclip_coverage.unwrap();

                let likelihood_max_softclip = |max_softclip| {
                    let mut lh = LogProb::ln_one();
                    for s in softclip_range() {
                        let count = common_obs.softclip_count(s);
                        let mu = if s <= max_softclip { varcov } else { 0.0 };
                        lh += poisson_pmf(count, mu);
                        if lh == LogProb::ln_zero() {
                            // stop early if we reach probability zero
                            break;
                        }
                    }
                    assert!(lh.is_valid());

                    lh
                };

                // calculate total probability to sample alt allele given the max_softclip
                // distribution
                let likelihoods: Vec<LogProb> = softclip_range().map(
                    |s| likelihood_max_softclip(s)
                ).collect_vec();
                let marginal = LogProb::ln_sum_exp(&likelihoods);
                assert!(!marginal.is_nan());
                assert!(marginal != LogProb::ln_zero(), "bug: marginal softclip dist prob of zero");

                let p = LogProb::ln_sum_exp(&probs.iter().enumerate().map(
                    |(max_softclip, prob_sample_alt)| {
                        // posterior probability for maximum softclip s
                        let prob_max_softclip = likelihoods[max_softclip as usize] - marginal;
                        // probability to sample from alt allele with this maximum softclip
                        prob_max_softclip + prob_sample_alt
                    }
                ).collect_vec()).cap_numerical_overshoot(NUMERICAL_EPSILON);
                assert!(p.is_valid(), "invalid probability {:?}", p);

                p
            }
        }
    }
}


/// An observation for or against a variant.
#[derive(Clone, Debug)]
pub struct Observation {
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    pub prob_mapping: LogProb,
    /// Probability that the read/read-pair comes from the alternative allele.
    pub prob_alt: LogProb,
    /// Probability that the read/read-pair comes from the reference allele.
    pub prob_ref: LogProb,
    /// Probability to sample the alt allele
    pub prob_sample_alt: LogProb,
    /// Type of evidence.
    pub evidence: Evidence
}


impl Observation {
    pub fn new(
        prob_mapping: LogProb,
        prob_alt: LogProb,
        prob_ref: LogProb,
        prob_sample_alt: LogProb,
        evidence: Evidence
    ) -> Self {
        Observation {
            prob_mapping: prob_mapping,
            prob_alt: prob_alt,
            prob_ref: prob_ref,
            prob_sample_alt: prob_sample_alt,
            evidence: evidence,
        }
    }

    pub fn is_alignment_evidence(&self) -> bool {
        if let Evidence::Alignment(_) = self.evidence {
            true
        } else {
            false
        }
    }
}


pub fn poisson_pmf(count: u32, mu: f64) -> LogProb {
    if mu == 0.0 {
        if count == 0 {
            LogProb::ln_one()
        } else {
            LogProb::ln_zero()
        }
    } else {
        LogProb(poisson_pdf(count, mu).ln())
    }
}


impl Serialize for Observation {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: Serializer {
        let mut s = serializer.serialize_struct("Observation", 3)?;
        s.serialize_field("prob_mapping", &self.prob_mapping)?;
        s.serialize_field("prob_alt", &self.prob_alt)?;
        s.serialize_field("prob_ref", &self.prob_ref)?;
        s.serialize_field("prob_sample_alt", &self.prob_sample_alt)?;
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


#[derive(Default, Clone, Debug)]
pub struct SoftclipObservation {
    counts: VecMap<u32>,
    start: Option<u32>,
    end: Option<u32>
}


impl SoftclipObservation {
    pub fn insert(&mut self, pos: u32, softclip: u32) {
        *self.counts.entry(softclip as usize).or_insert(0) += 1;
        self.start = Some(self.start.map_or(pos, |s| cmp::min(s, pos)));
        self.end = Some(self.end.map_or(pos, |e| cmp::max(e, pos)));
    }

    pub fn count(&self, softclip: u32) -> u32 {
        self.counts.get(softclip as usize).cloned().unwrap_or(0)
    }

    pub fn interval_len(&self) ->  u32 {
        self.end.map_or(0, |e| e - self.start.unwrap())
    }

    pub fn total_count(&self) -> u32 {
        self.counts.values().sum()
    }

    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }
}


/// Data that is shared among all observations over a locus.
#[derive(Clone, Debug)]
pub struct Common {
    pub leading_softclip_obs: Option<SoftclipObservation>,
    pub trailing_softclip_obs: Option<SoftclipObservation>,
    /// Average number of reads starting at any position in the region.
    pub softclip_coverage: Option<f64>,
    pub max_read_len: u32,
    pub enclosing_possible: bool
}


impl Common {
    pub fn new(variant: &Variant) -> Self {
        let enclosing_possible = variant.is_snv();

        Common {
            leading_softclip_obs: None,
            trailing_softclip_obs: None,
            max_read_len: 0,
            softclip_coverage: None,
            enclosing_possible: enclosing_possible
        }
    }

    pub fn softclip_count(&self, softclip: u32) -> u32 {
        if self.enclosing_possible {
            panic!("invalid operation: enclosing possible and softclips accessed");
        }
        self.leading_softclip_obs.as_ref().unwrap().count(softclip) +
        self.trailing_softclip_obs.as_ref().unwrap().count(softclip)
    }

    pub fn not_enough_softclips(&self) -> bool {
        if self.enclosing_possible {
            panic!("invalid operation: enclosing possible and softclips accessed");
        }
        self.leading_softclip_obs.as_ref().unwrap().total_count() < 2 ||
        self.trailing_softclip_obs.as_ref().unwrap().total_count() < 2
    }

    /// Update common observation with new reads.
    pub fn update(
        &mut self, records: &RecordBuffer, start: u32, variant: &Variant
    ) -> Result<(), Box<Error>> {
        let valid_records = || records.iter().filter(|rec| !rec.is_supplementary());

        // obtain maximum read len
        self.max_read_len = cmp::max(
            self.max_read_len,
            valid_records().map(
                |rec| rec.seq().len()
            ).max().unwrap_or(0) as u32
        );

        // determine if variant can be enclosed
        self.enclosing_possible |= {
            let max_indel_cigar = valid_records().map(|rec| {
                evidence::max_indel(&rec.cigar())
            }).max().unwrap_or(0);
            variant.len() <= max_indel_cigar
        };

        if !self.enclosing_possible {
            let leading_softclip_obs = self.leading_softclip_obs.get_or_insert_with(
                || SoftclipObservation::default()
            );
            let trailing_softclip_obs = self.trailing_softclip_obs.get_or_insert_with(
                || SoftclipObservation::default()
            );

            for rec in valid_records() {
                let cigar = rec.cigar();

                let leading_soft = evidence::Clips::leading(&cigar).soft();
                let trailing_soft = evidence::Clips::trailing(&cigar).soft();

                let overlap = Overlap::new(
                    rec, &cigar, start, variant, false, true
                )?;

                if overlap.is_none() {
                    continue;
                }
                if leading_soft > 0 || trailing_soft > 0 {
                    if leading_soft > trailing_soft {
                        leading_softclip_obs.insert(
                            (rec.pos() as u32).saturating_sub(leading_soft),
                            leading_soft
                        );
                    } else {
                        trailing_softclip_obs.insert(
                            (rec.pos() as u32).saturating_sub(leading_soft),
                            trailing_soft
                        );
                    }
                }
            }

            // update coverage
            self.softclip_coverage = if leading_softclip_obs.total_count() >= 2 ||
                                        trailing_softclip_obs.total_count() >= 2 {
                Some((
                    leading_softclip_obs.total_count() +
                    trailing_softclip_obs.total_count()
                ) as f64 / (
                    leading_softclip_obs.interval_len() + trailing_softclip_obs.interval_len()
                ) as f64)
            } else {
                None
            };

        } else {
            self.leading_softclip_obs = None;
            self.trailing_softclip_obs = None;
            self.softclip_coverage = None;
        }

        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_prob_sample_alt() {
    //     let obs = Observation {
    //         prob_mapping: LogProb::ln_one(),
    //         prob_alt: LogProb::ln_one(),
    //         prob_ref: LogProb::ln_one(),
    //         evidence: Evidence::dummy_alignment(),
    //         common: Common {
    //             softclip_obs: Some(softclip_obs),
    //             coverage:
    //         }
    //     }
    // }
}

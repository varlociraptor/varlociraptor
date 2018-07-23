use std::str;
use std::error::Error;
use std::cmp;
use std::ops::Range;
use std::f64;
use std::collections::BTreeMap;

use serde::Serialize;
use serde::ser::{Serializer, SerializeStruct};
use itertools::Itertools;
use rgsl::randist::poisson::poisson_pdf;

use bio::stats::LogProb;
use rust_htslib::bam;
use rust_htslib::bam::record::{CigarString, Cigar};

use model::Variant;
use model::evidence;
use model::sample::RecordBuffer;
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
    /// Otherwise, this has three components.
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
    /// softclips that overlap the variant and dividing by the interval length.
    /// Then, we can calculate the likelihood of the observed softclips given an assumed maximum
    /// softclip by taking the product of the poisson distributed probabilities for each observed
    /// softclip count.
    ///
    /// # Component 2: feasible indel cigar positions
    /// Indel operations in cigar strings tend to be harder to obtain towards the end of the read.
    /// Hence, difficult (longer) indels are often only observed towards the center of the read,
    /// leading to an underestimation of allele frequencies. We again assume poisson distributed
    /// reads, calculate the mean over a given assumed feasible interval and obtain the likelihood
    /// for a given feasible interval by a product of poisson distributed probabilities.
    /// Then, we calculate the overall likelihood for a given number of feasible positions considering
    /// both softclips and indel positions via dynamic programming.
    /// By applying Bayes theorem, we obtain a posterior probability for each number of feasible
    /// positions.
    ///
    /// # Component 3:
    /// We calculate the probability to sample a fragment from the alt allele given a number of
    /// feasible positions and the read lengths.
    ///
    /// # Joint probability
    /// The final result is obtained by combining the two components to a total probability.
    ///
    /// # Arguments
    /// * allele_freq - Given allele frequency of sample.
    pub fn joint_prob(&self, common_obs: &Common) -> LogProb {
        match self {
            &ProbSampleAlt::One => LogProb::ln_one(),
            &ProbSampleAlt::Independent(p) => p,
            &ProbSampleAlt::Dependent(ref probs) => {
                let prob_feasible = common_obs.prob_feasible.as_ref().unwrap();

                let p = LogProb::ln_sum_exp(&probs.iter().enumerate().map(
                    |(feasible, prob_sample_alt)| {
                        // probability to sample from alt allele with this number of feasible
                        // positions
                        prob_feasible[feasible as usize] + prob_sample_alt
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


type CigarObservations = BTreeMap<u32, u32>;


/// Data that is shared among all observations over a locus.
#[derive(Clone, Debug, Default)]
pub struct Common {
    softclip_obs: BTreeMap<u32, u32>,
    indel_obs: BTreeMap<u32, u32>,
    pub max_read_len: u32,
    prob_feasible: Option<Vec<LogProb>>
}


impl Common {

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

        match variant {
            &Variant::Deletion(_) | &Variant::Insertion(_) => {
                for rec in valid_records() {
                    let cigar = rec.cigar_cached().unwrap();

                    let overlap = Overlap::new(
                        rec, cigar, start, variant, true
                    )?;

                    if overlap.is_none() {
                        continue;
                    }

                    // record softclips
                    let leading_soft = evidence::Clips::leading(&cigar).soft();
                    let trailing_soft = evidence::Clips::trailing(&cigar).soft();
                    if leading_soft > 0 {
                        *self.softclip_obs.entry(leading_soft).or_insert(0) += 1;
                    }
                    if trailing_soft > 0 {
                        *self.softclip_obs.entry(trailing_soft).or_insert(0) += 1;
                    }

                    // record indel operations
                    let mut qpos = 0;
                    for c in cigar {
                        match c {
                            &Cigar::Del(l) => {
                                if let &Variant::Deletion(m) = variant {
                                    if l >= m {
                                        *self.indel_obs.entry(qpos).or_insert(0) += 1;
                                    }
                                }
                                // do not increase qpos in case of a deletion
                            },
                            &Cigar::Ins(l) => {
                                if let &Variant::Insertion(ref seq) = variant {
                                    if l >= seq.len() as u32 {
                                        *self.indel_obs.entry(qpos).or_insert(0) += 1;
                                    }
                                }
                                qpos += l;
                            },
                            &Cigar::SoftClip(l) | &Cigar::Match(l) | &Cigar::Equal(l) |
                            &Cigar::Diff(l) | &Cigar::Pad(l) => {
                                qpos += l;
                            },
                            &Cigar::RefSkip(_) | &Cigar::HardClip(_) => {
                                // do nothing
                            }
                        }

                        assert!(
                            qpos <= rec.seq().len() as u32,
                            format!("bug: qpos larger than read len ({})", cigar)
                        );
                    }
                }
            },
            &Variant::SNV(_) | &Variant::None => ()
        }

        Ok(())
    }

    pub fn finalize(&mut self) {
        if self.softclip_obs.is_empty() && self.indel_obs.is_empty() {
            let mut prob_feasible = vec![LogProb::ln_zero(); self.max_read_len as usize + 1];
            prob_feasible[0] = LogProb::ln_one();
            self.prob_feasible = Some(prob_feasible);
            return;
        }

        let likelihoods = {

            let likelihood_softclip = if !self.softclip_obs.is_empty() {
                let mut lhs = (0..self.max_read_len + 1).map(|max_softclip| {
                    self.likelihood_feasible(0..max_softclip, &self.softclip_obs)
                }).collect_vec();

                if lhs.iter().all(|lh| *lh == LogProb::ln_zero()) {
                    // There is apparently no max softclip that can explain the observed counts
                    // when assuming they are poisson distributed.
                    // The best we can do now is to simply take the maximum observation.
                    let max_softclip = *self.softclip_obs.keys().max().unwrap();
                    lhs[max_softclip as usize] = LogProb::ln_one();
                }

                lhs
            } else {
                let mut l = vec![LogProb::ln_zero(); self.max_read_len as usize + 1];
                l[0] = LogProb::ln_one();
                l
            };

            if !self.indel_obs.is_empty() {
                let mut likelihoods = vec![LogProb::ln_zero(); self.max_read_len as usize + 1];
                for start in 0..self.max_read_len / 2 {
                    let end = self.max_read_len - start + 1;
                    let lh = self.likelihood_feasible(start..end, &self.indel_obs);
                    for max_softclip in 0..self.max_read_len + 1 {
                        let infeasible = start.saturating_sub(max_softclip) + (self.max_read_len - (end - 1));
                        let f = (self.max_read_len as usize).saturating_sub(infeasible as usize);
                        likelihoods[f] = likelihoods[f].ln_add_exp(
                            lh + likelihood_softclip[max_softclip as usize]
                        );
                    }
                }
                likelihoods
            } else {
                likelihood_softclip
            }
        };

        // calculate joint probability to sample alt allele given the feasible position
        // distribution
        let marginal = LogProb::ln_sum_exp(&likelihoods);
        assert!(!marginal.is_nan());
        assert!(marginal != LogProb::ln_zero(), "bug: marginal feasibility dist prob of zero");

        self.prob_feasible = Some(likelihoods.iter().map(|lh| lh - marginal).collect_vec());
    }

    fn likelihood_feasible(
        &self, range: Range<u32>, obs: &CigarObservations
    ) -> LogProb {
        assert!(range.start <= range.end);
        let coverage = if range.start == range.end {
            0.0
        } else {
            obs.range(range.clone()).map(|(_, count)| count).sum::<u32>() as f64 /
            (range.end - range.start) as f64
        };
        assert!(coverage != f64::INFINITY);

        let mut lh = LogProb::ln_one();
        for s in 0..self.max_read_len as u32 + 1 {
            let count = *obs.get(&s).unwrap_or(&0);
            let mu = if s >= range.start && s < range.end { coverage } else { 0.0 };
            lh += poisson_pmf(count, mu);
            if lh == LogProb::ln_zero() {
                // stop early if we reach probability zero
                break;
            }
        }
        assert!(lh.is_valid());

        lh
    }
}

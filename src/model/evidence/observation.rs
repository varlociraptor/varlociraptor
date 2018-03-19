use std::str;
use std::error::Error;
use std::cmp;
use std::ops::Range;

use vec_map::VecMap;
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

                if common_obs.total_indel_count() == 1 {
                    // we have only one indel. Since we cannot infer a distribution from this,
                    // we assume that all positions are feasible.
                    return *probs.iter().last().unwrap();
                }

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


#[derive(Default, Clone, Debug)]
pub struct CigarObservation {
    counts: VecMap<u32>,
    start: Option<u32>,
    end: Option<u32>
}


impl CigarObservation {
    pub fn insert(&mut self, pos: u32, observation: u32) {
        *self.counts.entry(observation as usize).or_insert(0) += 1;
        self.start = Some(self.start.map_or(pos, |s| cmp::min(s, pos)));
        self.end = Some(self.end.map_or(pos, |e| cmp::max(e, pos)));
    }

    pub fn count(&self, observation: u32) -> u32 {
        self.counts.get(observation as usize).cloned().unwrap_or(0)
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

    /// Maximum observation, zero if nothing observed.
    pub fn max_observation(&self) -> Option<u32> {
        self.counts.keys().last().map(|o| o as u32)
    }

    /// Minimum observation, zero if nothing observed.
    pub fn min_observation(&self) -> Option<u32> {
        self.counts.keys().next().map(|o| o as u32)
    }
}


/// Data that is shared among all observations over a locus.
#[derive(Clone, Debug, Default)]
pub struct Common {
    pub leading_softclip_obs: Option<CigarObservation>,
    pub trailing_softclip_obs: Option<CigarObservation>,
    pub indel_obs: Option<CigarObservation>,
    /// Average number of reads starting at any position in the region.
    pub softclip_coverage: Option<f64>,
    pub indel_coverage: Option<f64>,
    pub max_read_len: u32,
    prob_feasible: Option<Vec<LogProb>>
}


impl Common {
    pub fn indel_count(&self, qpos: u32) -> u32 {
        self.indel_obs.as_ref().unwrap().count(qpos)
    }

    pub fn softclip_count(&self, softclip: u32) -> u32 {
        self.leading_softclip_obs.as_ref().unwrap().count(softclip) +
        self.trailing_softclip_obs.as_ref().unwrap().count(softclip)
    }

    pub fn not_enough_softclips(&self) -> bool {
        self.leading_softclip_obs.as_ref().unwrap().total_count() < 2 ||
        self.trailing_softclip_obs.as_ref().unwrap().total_count() < 2
    }

    pub fn total_indel_count(&self) -> u32 {
        self.indel_obs.as_ref().unwrap().total_count()
    }

    pub fn max_softclip(&self) -> u32 {
        cmp::max(
            self.leading_softclip_obs.as_ref().unwrap().max_observation().unwrap_or(0),
            self.trailing_softclip_obs.as_ref().unwrap().max_observation().unwrap_or(0)
        )
    }

    pub fn min_indel_pos(&self) -> u32 {
        self.indel_obs.as_ref().unwrap().start.unwrap()
    }

    pub fn max_indel_pos(&self) -> u32 {
        self.indel_obs.as_ref().unwrap().end.unwrap()
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

        match variant {
            &Variant::Deletion(_) | &Variant::Insertion(_) => {
                let leading_softclip_obs = self.leading_softclip_obs.get_or_insert_with(
                    || CigarObservation::default()
                );
                let trailing_softclip_obs = self.trailing_softclip_obs.get_or_insert_with(
                    || CigarObservation::default()
                );
                let indel_obs = self.indel_obs.get_or_insert_with(
                    || CigarObservation::default()
                );

                for rec in valid_records() {
                    let cigar = rec.cigar();

                    let overlap = Overlap::new(
                        rec, &cigar, start, variant, true
                    )?;

                    if overlap.is_none() {
                        continue;
                    }

                    // record softclips
                    let leading_soft = evidence::Clips::leading(&cigar).soft();
                    let trailing_soft = evidence::Clips::trailing(&cigar).soft();
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

                    // record indel operations
                    let mut qpos = 0;
                    for c in &cigar {
                        match (c, variant) {
                            (&Cigar::Del(l), &Variant::Deletion(m)) if l >= m => {
                                indel_obs.insert(qpos, qpos);
                            },
                            (&Cigar::Ins(l), &Variant::Insertion(ref seq)) if l >= seq.len() as u32 => {
                                indel_obs.insert(qpos, qpos);
                            },
                            _ => ()
                        }

                        qpos += c.len();
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

                self.indel_coverage = if indel_obs.total_count() >= 2 {
                    Some(indel_obs.total_count() as f64 / indel_obs.interval_len() as f64)
                } else {
                    None
                };
            },
            &Variant::SNV(_) | &Variant::None => {
                self.leading_softclip_obs = None;
                self.trailing_softclip_obs = None;
                self.softclip_coverage = None;
                self.indel_coverage = None;
            }
        }

        Ok(())
    }

    pub fn finalize(&mut self) {
        if self.softclip_coverage.is_none() && self.indel_coverage.is_none() {
            let mut prob_feasible = vec![LogProb::ln_zero(); self.max_read_len as usize + 1];
            prob_feasible[0] = LogProb::ln_one();
            self.prob_feasible = Some(prob_feasible);
            return;
        }

        let likelihoods = {
            let feasible_range = || 0..self.max_read_len as u32 + 1;

            let likelihood_softclip = if !self.not_enough_softclips() {
                let softclip_cov = self.softclip_coverage.unwrap();
                feasible_range().map(|max_softclip| {
                    Self::likelihood_feasible(|s| s <= max_softclip, softclip_cov, |s| self.softclip_count(s), feasible_range())
                }).collect_vec()
            } else {
                // Not enough softclips observed.
                // We take the maximum observed softclip as the true maximum.
                let s = self.max_softclip() as usize;
                let mut likelihoods = vec![LogProb::ln_zero(); self.max_read_len as usize + 1];
                likelihoods[s] = LogProb::ln_one();
                likelihoods
            };

            if let Some(indel_cov) = self.indel_coverage {
                let mut likelihoods = vec![LogProb::ln_zero(); self.max_read_len as usize + 1];
                for start in 0..self.min_indel_pos() {
                    for end in cmp::max(start, self.max_indel_pos())..self.max_read_len as u32 {
                        let lh = Self::likelihood_feasible(
                            |pos| pos >= start && pos < end,
                            indel_cov,
                            |pos| self.indel_count(pos),
                            feasible_range()
                        );
                        for max_softclip in feasible_range() {
                            let infeasible = start.saturating_sub(
                                max_softclip
                            ) + (self.max_read_len as u32).saturating_sub(
                                max_softclip
                            ).saturating_sub(end);
                            let f = (self.max_read_len as usize).saturating_sub(infeasible as usize);
                            likelihoods[f] = likelihoods[f].ln_add_exp(
                                lh + likelihood_softclip[max_softclip as usize]
                            );
                        }
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
        assert!(marginal != LogProb::ln_zero(), "bug: marginal softclip dist prob of zero");

        self.prob_feasible = Some(likelihoods.iter().map(|lh| lh - marginal).collect_vec());
    }

    fn likelihood_feasible<V: Fn(u32) -> bool, C: Fn(u32) -> u32>(
        is_valid: V, varcov: f64, count: C, feasible_range: Range<u32>
    ) -> LogProb {
        let mut lh = LogProb::ln_one();
        for s in feasible_range {
            //let count = common_obs.softclip_count(s);
            let count = count(s);
            //let mu = if s <= max_softclip { varcov } else { 0.0 };
            let mu = if is_valid(s) { varcov } else { 0.0 };
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

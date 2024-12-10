// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::ops::{AddAssign, Deref};
use std::path::PathBuf;
use std::str;

use anyhow::Result;
use bio::stats::{
    bayesian::bayes_factors::evidence::KassRaftery, bayesian::bayes_factors::BayesFactor, LogProb,
    PHREDProb, Prob,
};
use counter::Counter;
use half::f16;
use itertools::join;
use itertools::Itertools;
use ordered_float::NotNan;
use rust_htslib::bcf::Read;
use rust_htslib::{bam, bam::record::Cigar, bcf};
use serde::Serializer;

use crate::variants::model;
use crate::Event;

pub(crate) mod adaptive_integration;
pub(crate) mod anonymize;
pub(crate) mod aux_info;
pub mod collect_variants;
pub(crate) mod comparison;
pub(crate) mod homopolymers;
pub(crate) mod log2_fold_change;
pub(crate) mod variant_buffer;

pub(crate) use collect_variants::collect_variants;

pub(crate) const NUMERICAL_EPSILON: f64 = 1e-3;

pub(crate) type PathMap = HashMap<String, PathBuf>;

lazy_static! {
    pub(crate) static ref PROB_05: LogProb = LogProb::from(Prob(0.5f64));
    pub(crate) static ref PROB_033: LogProb = LogProb::from(Prob(1.0 / 3.0));
    pub(crate) static ref PROB_025: LogProb = LogProb::from(Prob(0.25));
    pub(crate) static ref PROB_095: LogProb = LogProb::from(Prob(0.95));
    pub(crate) static ref PROB_09: LogProb = LogProb::from(Prob(0.9));
}

pub(crate) fn aux_tag_strand_info(record: &bam::Record) -> Option<&[u8]> {
    if let Ok(bam::record::Aux::String(strand_info)) = record.aux(b"SI") {
        Some(strand_info.as_bytes())
    } else {
        None
    }
}

pub(crate) fn aux_tag_is_entire_fragment(record: &bam::Record) -> bool {
    match record.aux(b"EF") {
        Ok(bam::record::Aux::U8(entire_fragment)) => entire_fragment == 1,
        Ok(bam::record::Aux::I8(entire_fragment)) => entire_fragment == 1,
        Ok(bam::record::Aux::U16(entire_fragment)) => entire_fragment == 1,
        Ok(bam::record::Aux::I16(entire_fragment)) => entire_fragment == 1,
        Ok(bam::record::Aux::U32(entire_fragment)) => entire_fragment == 1,
        Ok(bam::record::Aux::I32(entire_fragment)) => entire_fragment == 1,
        _ => false,
    }
}

/// Checks whether the given BCF contains fields required for evaluating haplotypes.
/// Currently, this means that the EVENT or the MATEID field has to be defined in the
/// header.
pub(crate) fn is_haplotype_bcf(reader: &bcf::Reader) -> bool {
    for rec in reader.header().header_records() {
        if let bcf::header::HeaderRecord::Info { values, .. } = rec {
            if values
                .get("ID")
                .map_or(false, |id| id == "EVENT" || id == "MATEID")
            {
                return true;
            }
        }
    }
    false
}

pub fn is_bnd(record: &mut bcf::Record) -> Result<bool> {
    Ok(record
        .info(b"SVTYPE")
        .string()?
        .map_or(false, |entries| entries[0] == b"BND"))
}

pub(crate) fn info_tag_svtype(record: &mut bcf::Record) -> Result<Option<Vec<u8>>> {
    Ok(record.info(b"SVTYPE").string()?.map(|v| v[0].to_owned()))
}

pub(crate) fn info_tag_event(record: &mut bcf::Record) -> Result<Option<Vec<u8>>> {
    Ok(record.info(b"EVENT").string()?.map(|v| v[0].to_owned()))
}

pub(crate) fn info_tag_mateid(record: &mut bcf::Record) -> Result<Option<Vec<u8>>> {
    // TODO support multiple mateids (in case of uncertainty, see spec)
    Ok(record.info(b"MATEID").string()?.map(|v| v[0].to_owned()))
}

pub(crate) fn is_reverse_strand(record: &bam::Record) -> bool {
    record.flags() & 0x10 != 0
}

pub(crate) fn contains_indel_op(record: &bam::Record) -> bool {
    record
        .cigar_cached()
        .expect("bug: cigar accessed before caching")
        .iter()
        .any(|op| matches!(op, Cigar::Ins(_) | Cigar::Del(_)))
}

pub(crate) fn generalized_cigar<T: Hash + Eq + Clone + Display, F, K>(
    items: impl Iterator<Item = T>,
    keep_order: bool,
    aux_sort: F,
) -> String
where
    F: FnMut(&(T, usize)) -> K,
    K: Ord,
{
    if keep_order {
        join(
            items
                .map(|item| (item, 1))
                .coalesce(|(a, n), (b, m)| {
                    if a == b {
                        Ok((a, n + m))
                    } else {
                        Err(((a, n), (b, m)))
                    }
                })
                .map(|(item, count)| format!("{}{}", count, item)),
            "",
        )
    } else {
        let items: Counter<T> = items.collect();
        join(
            items
                .most_common()
                .into_iter()
                .sorted_by_key(aux_sort)
                .map(|(item, count)| format!("{}{}", count, item)),
            "",
        )
    }
}

pub(crate) fn bayes_factor_to_letter(bayes_factor: BayesFactor) -> char {
    match bayes_factor.evidence_kass_raftery() {
        KassRaftery::Barely => 'B',
        KassRaftery::None if relative_eq!(*bayes_factor, 1.0) => 'E',
        KassRaftery::None => 'N',
        KassRaftery::Positive => 'P',
        KassRaftery::Strong => 'S',
        KassRaftery::VeryStrong => 'V',
    }
}

/// Sum up in log space the probabilities of the given tags for all variants of
/// vartype in the given BCF record.
///
/// # Arguments
///
/// * `record` - BCF record
/// * `tags` - tags of the set of events to sum up for a particular site and variant
/// * `vartype` - the variant type to consider
pub(crate) fn tags_prob_sum(
    record: &mut bcf::Record,
    tags: &[String],
    vartype: Option<&model::VariantType>,
) -> Result<Vec<Option<LogProb>>> {
    let mut skips = SimpleCounter::default();
    let variants = collect_variants(record, false, Some(&mut skips))?;
    let mut tags_probs_out = vec![Vec::new(); variants.len()];

    for tag in tags {
        if let Some(tags_probs_in) = (record.info(tag.as_bytes()).float())? {
            //tag present
            for (i, (variant_info, tag_prob)) in
                variants.iter().zip(tags_probs_in.iter()).enumerate()
            {
                if (vartype.is_some() && !variant_info.variant().is_type(vartype.unwrap()))
                    || tag_prob.is_nan()
                {
                    continue;
                }
                tags_probs_out[i].push(LogProb::from(PHREDProb(*tag_prob as f64)));
            }
        }
    }

    Ok(tags_probs_out
        .into_iter()
        .map(|probs| {
            if !probs.is_empty() {
                Some(LogProb::ln_sum_exp(&probs).cap_numerical_overshoot(NUMERICAL_EPSILON))
            } else {
                None
            }
        })
        .collect_vec())
}

pub(crate) fn events_to_tags<E>(events: &[E]) -> Vec<String>
where
    E: Event,
{
    events.iter().map(|e| event_to_tag(e)).collect_vec()
}

pub(crate) fn event_to_tag<E>(event: &E) -> String
where
    E: Event,
{
    event.tag_name("PROB")
}

/// Collect distribution of posterior probabilities from a VCF file that has been written by
/// varlociraptor.
///
/// # Arguments
///
/// * `calls` - BCF reader with varlociraptor calls
/// * `events` - the set of events to sum up for a particular site
/// * `vartype` - the variant type to consider
pub(crate) fn collect_prob_dist<E>(
    calls: &mut bcf::Reader,
    events: &[E],
    vartype: Option<&model::VariantType>,
) -> Result<Vec<NotNan<f64>>>
where
    E: Event,
{
    let mut visited_breakend_events = HashSet::new();

    let mut record = calls.empty_record();
    let mut prob_dist = Vec::new();
    let tags = events_to_tags(events);
    loop {
        match calls.read(&mut record) {
            None => break,
            Some(res) => res?,
        }
        if let Ok(Some(event)) = info_tag_event(&mut record) {
            if visited_breakend_events.contains(&event) {
                // Do not record probability of this event twice.
                continue;
            } else {
                visited_breakend_events.insert(event.to_owned());
            }
        }

        for p in (tags_prob_sum(&mut record, &tags, vartype)?)
            .into_iter()
            .flatten()
        {
            prob_dist.push(NotNan::new(*p)?);
        }
    }
    prob_dist.sort();
    Ok(prob_dist)
}

/// Filter a VCF record stream by a minimum threshold on the sum of
/// posterior probabilities of a given set of Events. The threshold
/// should be an informative false discovery rate (FDR) threshold,
/// e.g. determined with the varlociraptor FDR control functionality.
///
/// # Arguments
///
/// * `calls` - BCF reader with varlociraptor calls
/// * `threshold` - minimum threshold for the sum of posterior probabilities of the set of Events considered
/// * `calls` - BCF writer for the filtered varlociraptor calls
/// * `events` - the set of Events to filter on
/// * `vartype` - the variant type to consider (if None, use all types)
/// * `smart` - if true, use a smart thresholding strategy that just filters based on presence with the threshold
///             and then keeps records whether the sum of the event probabilities exceeds 50%.
pub(crate) fn filter_by_threshold<E: Event>(
    calls: &mut bcf::Reader,
    threshold: Option<LogProb>,
    out: &mut bcf::Writer,
    events: &[E],
    vartype: Option<&model::VariantType>,
    smart: bool,
) -> Result<()> {
    let mut breakend_event_decisions = HashMap::new();

    let tags = events_to_tags(events);
    let absent_and_artifact_tags = ["PROB_ABSENT".to_owned(), "PROB_ARTIFACT".to_owned()];
    let filter = |record: &mut bcf::Record| -> Result<Vec<bool>> {
        let bnd_event = info_tag_event(record).ok().flatten();
        let keep = if let Some(event) = bnd_event.as_ref() {
            breakend_event_decisions.get(event).cloned()
        } else {
            None
        };

        let probs_events = tags_prob_sum(record, &tags, vartype)?;
        let probs_absent_or_artifact = if smart {
            tags_prob_sum(record, &absent_and_artifact_tags, vartype)?
        } else {
            vec![None; probs_events.len()]
        };

        assert!(
            bnd_event.is_none() || probs_events.len() == 1,
            "breakend events may only occur in single variant records"
        );

        Ok(probs_events
            .into_iter()
            .zip(probs_absent_or_artifact)
            .map(|(prob_events, prob_absent_or_artifact)| {
                if let Some(keep) = keep {
                    // already know decision from previous breakend
                    keep
                } else {
                    let p = if smart {
                        // METHOD: in smart mode, the record is kept if the probability for presence exceeds the threshold
                        prob_absent_or_artifact.map(|p| p.ln_one_minus_exp())
                    } else {
                        prob_events
                    };

                    let mut keep = match (p, threshold) {
                        // we allow some numerical instability in case of equality
                        (Some(p), Some(threshold))
                            if p > threshold || relative_eq!(*p, *threshold) =>
                        {
                            true
                        }
                        (Some(_), None) => true,
                        _ => false,
                    };

                    if smart {
                        // METHOD: in smart mode, we only keep the record if the sum of the event probabilities
                        // exceeds 50%.
                        keep &= prob_events.map(|p| p > *PROB_05).unwrap_or(false);
                    }

                    if let Some(event) = bnd_event.as_ref() {
                        breakend_event_decisions.insert(event.to_owned(), keep);
                    }

                    keep
                }
            })
            .collect())
    };
    filter_calls(calls, out, filter)
}

/// Filter calls by a given function.
///
/// # Arguments
/// * `calls` - the calls to filter
/// * `out` - output BCF
/// * `filter` - function to filter by. Has to return a bool for every alternative allele.
///   True means to keep the allele.
pub(crate) fn filter_calls<F, I, II>(
    calls: &mut bcf::Reader,
    out: &mut bcf::Writer,
    mut filter: F,
) -> Result<()>
where
    F: FnMut(&mut bcf::Record) -> Result<II>,
    I: Iterator<Item = bool>,
    II: IntoIterator<Item = bool, IntoIter = I>,
{
    let mut record = calls.empty_record();
    let mut i = 1;
    loop {
        match calls.read(&mut record) {
            None => return Ok(()),
            Some(res) => res?,
        }

        let mut remove = vec![false]; // don't remove the reference allele
        remove.extend(filter(&mut record)?.into_iter().map(|keep| !keep));

        assert_eq!(
            remove.len(),
            record.allele_count() as usize,
            "bug: filter passed to filter_calls has to return a bool for each alt allele at record {}.",
            i,
        );

        // Write trimmed record if any allele remains. Otherwise skip the record.
        if !remove[1..].iter().all(|r| *r) {
            record.remove_alleles(&remove)?;
            out.write(&record)?;
        }

        i += 1;
    }
}

/// Returns true if all PROB_{event}s are PHRED scaled
pub(crate) fn is_phred_scaled(inbcf: &bcf::Reader) -> bool {
    get_event_tags(inbcf)
        .iter()
        // check for missing closing parenthesis for backward compatibility
        .all(|(_, d)| d.ends_with("(PHRED)") || !d.ends_with(')'))
}

/// Returns (ID, Description) for each PROB_{event} INFO tag
pub(crate) fn get_event_tags(inbcf: &bcf::Reader) -> Vec<(String, String)> {
    inbcf
        .header()
        .header_records()
        .into_iter()
        .filter_map(|rec| {
            if let bcf::header::HeaderRecord::Info { values, .. } = rec {
                let id = values["ID"].clone();
                if id.starts_with("PROB_") {
                    let description = values["Description"].clone();
                    let description = description.trim_matches('"').into();
                    return Some((id, description));
                }
            }
            None
        })
        .collect_vec()
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub(crate) enum MiniLogProb {
    F16(f16),
    F32(f32),
}

impl MiniLogProb {
    /// Convert LogProb into a minimal representation for storage.
    /// If integer part is less than -1 and can be represented in f16,
    /// we use f16. Else, we use f32.
    pub(crate) fn new(prob: LogProb) -> Self {
        let half = f16::from_f64(*prob);
        let proj = half.to_f64();
        if *prob < -10.0 && proj.floor() as i64 == prob.floor() as i64 {
            MiniLogProb::F16(half)
        } else {
            MiniLogProb::F32(*prob as f32)
        }
    }

    pub(crate) fn to_logprob(&self) -> LogProb {
        LogProb(match self {
            MiniLogProb::F16(p) => p.to_f64(),
            MiniLogProb::F32(p) => *p as f64,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::variants::model::VariantType;
    use crate::SimpleEvent;
    use bio::stats::{LogProb, Prob};
    use rust_htslib::bcf::{self, Read};

    #[test]
    fn test_tags_prob_sum() {
        // set up test input
        let test_file = "tests/resources/test_tags_prob_sum/overshoot.vcf";
        let mut overshoot_calls = bcf::Reader::from_path(test_file).unwrap();
        let mut record = overshoot_calls.empty_record();
        match overshoot_calls.read(&mut record) {
            None => panic!("BCF reading error: seems empty"),
            Some(res) => res.unwrap(),
        }

        // set up all alt events with names as in prosolo
        let alt_tags = [
            String::from("PROB_ADO_TO_REF"),
            String::from("PROB_ADO_TO_ALT"),
            String::from("PROB_HOM_ALT"),
            String::from("PROB_HET"),
            String::from("PROB_ERR_REF"),
        ];

        let snv = VariantType::Snv;

        if let Ok(prob_sum) = tags_prob_sum(&mut record, &alt_tags, Some(&snv)) {
            assert_eq!(LogProb::ln_one(), prob_sum[0].unwrap());
        } else {
            panic!("tags_prob_sum(&overshoot_calls, &alt_events, &snv) returned Error")
        }
    }

    #[test]
    fn test_collect_prob_dist() {
        // setup events with names as in varlociraptor2
        let events = vec![
            SimpleEvent {
                name: "germline".to_owned(),
            },
            SimpleEvent {
                name: "somatic".to_owned(),
            },
        ];
        // setup absent event as the complement of the other events
        let absent_event = vec![SimpleEvent {
            name: "absent".to_owned(),
        }];

        let test_file = "tests/resources/test_collect_prob_dist/min.calls.vcf";

        //TESTS deletion
        let del = VariantType::Deletion(None);

        let mut del_calls_1 = bcf::Reader::from_path(test_file).unwrap();
        let prob_del = collect_prob_dist(&mut del_calls_1, &events, Some(&del)).unwrap();
        println!("prob_del[0]: {:?}", prob_del[0].into_inner());
        assert_eq!(prob_del.len(), 1);
        assert_relative_eq!(prob_del[0].into_inner(), Prob(0.8).ln(), epsilon = 0.000005);

        let mut del_calls_2 = bcf::Reader::from_path(test_file).unwrap();
        let prob_del_abs = collect_prob_dist(&mut del_calls_2, &absent_event, Some(&del)).unwrap();
        assert_eq!(prob_del_abs.len(), 1);
        assert_relative_eq!(
            prob_del_abs[0].into_inner(),
            Prob(0.2).ln(),
            epsilon = 0.000005
        );

        //TESTS insertion
        let ins = VariantType::Insertion(None);

        let mut ins_calls_1 = bcf::Reader::from_path(test_file).unwrap();
        let prob_ins = collect_prob_dist(&mut ins_calls_1, &events, Some(&ins)).unwrap();
        assert_eq!(prob_ins.len(), 1);
        assert_relative_eq!(prob_ins[0].into_inner(), Prob(0.2).ln(), epsilon = 0.000005);

        let mut ins_calls_2 = bcf::Reader::from_path(test_file).unwrap();
        let prob_ins_abs = collect_prob_dist(&mut ins_calls_2, &absent_event, Some(&ins)).unwrap();
        assert_eq!(prob_ins_abs.len(), 1);
        assert_relative_eq!(
            prob_ins_abs[0].into_inner(),
            Prob(0.8).ln(),
            epsilon = 0.000005
        );
    }

    #[test]
    fn test_filter_by_threshold() {
        // TODO: make this test work with both thresholds, testing against expected_output files
        /*
        // set up test input
        let test_file = "tests/resources/test_tags_prob_sum/overshoot.vcf";
        let mut calls = bcf::Reader::from_path( test_file ).unwrap();

        let threshold_1 = 0.1;
        let threshold_2 = 0.00000000001;

        let events = vec![
            SimpleEvent { name: "ADO_TO_REF".to_owned() },
            SimpleEvent { name: "ADO_TO_ALT".to_owned() },
            SimpleEvent { name: "HOM_ALT".to_owned() },
            SimpleEvent { name: "HET".to_owned() },
            SimpleEvent { name: "ERR_REF".to_owned() }
        ];

        let snv = VariantType::SNV;

        let header = bcf::Header::with_template(&calls.header());
        let mut out = bcf::Writer::from_stdout(&header, false, false).unwrap();

        filter_by_threshold(&mut calls, &threshold, &mut out, &events, &snv);

        panic!("Just checking");
        */
    }
}

#[derive(CopyGetters, Debug, Serialize, Deserialize, Clone)]
pub struct SimpleCounter<T>
where
    T: Eq + Hash + Debug,
{
    #[serde(serialize_with = "serialize_tuple_keys")]
    inner: HashMap<T, usize>,
    #[getset(get_copy = "pub(crate)")]
    total_count: usize,
}

fn serialize_tuple_keys<T: Eq + Hash + Debug, S>(
    map: &HashMap<T, usize>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.collect_map(map.iter().map(|(k, v)| (format!("{:?}", &k), v)))
}

impl<T> SimpleCounter<T>
where
    T: Eq + Hash + Debug,
{
    pub(crate) fn incr(&mut self, event: T) {
        self.total_count += 1;
        *self.inner.entry(event).or_insert(0) += 1;
    }
}

impl<T> Default for SimpleCounter<T>
where
    T: Eq + Hash + Debug,
{
    fn default() -> Self {
        SimpleCounter {
            inner: HashMap::new(),
            total_count: 0,
        }
    }
}

impl<T> Deref for SimpleCounter<T>
where
    T: Eq + Hash + Debug,
{
    type Target = HashMap<T, usize>;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T> AddAssign for SimpleCounter<T>
where
    T: Eq + Hash + Debug,
{
    fn add_assign(&mut self, rhs: Self) {
        self.total_count += rhs.total_count;
        rhs.inner.into_iter().for_each(move |(k, v)| {
            *self.inner.entry(k).or_insert(0) += v;
        });
    }
}

/// Format to given precision
/// Taken from https://stackoverflow.com/a/60498499/7070491
pub fn signif(float: f64, precision: usize) -> String {
    // compute absolute value
    let a = float.abs();

    // if abs value is greater than 1, then precision becomes less than "standard"
    let precision = if a >= 1. {
        // reduce by number of digits, minimum 0
        let n = (1. + a.log10().floor()) as usize;
        if n <= precision {
            precision - n
        } else {
            0
        }
    // if precision is less than 1 (but non-zero), then precision becomes greater than "standard"
    } else if a > 0. {
        // increase number of digits
        let n = -(1. + a.log10().floor()) as usize;
        precision + n
    // special case for 0
    } else {
        0
    };

    // format with the given computed precision
    format!("{0:.1$}", float, precision)
}

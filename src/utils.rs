// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::fmt::Display;
use std::fs;
use std::hash::Hash;
use std::ops::Range;
use std::str;

use anyhow::Result;
use bio::io::fasta;
use bio::stats::{bayesian::bayes_factors::evidence::KassRaftery, LogProb, PHREDProb};
use counter::Counter;
use half::f16;
use itertools::join;
use itertools::Itertools;
use ordered_float::NotNan;
use rust_htslib::bcf::Read;
use rust_htslib::{bam, bcf, bcf::record::Numeric};

use crate::errors;
use crate::model;
use crate::utils;
use crate::Event;

pub const NUMERICAL_EPSILON: f64 = 1e-3;

#[derive(new, Getters, CopyGetters, Debug)]
pub struct GenomicLocus {
    #[getset(get = "pub")]
    chrom: Vec<u8>,
    #[getset(get_copy = "pub")]
    pos: u32,
}

/// Select values with given indices from a slice and return them as an iterator.
pub fn select<'a, T: Clone>(idx: &'a [usize], values: &'a [T]) -> impl Iterator<Item = T> + 'a {
    idx.iter().map(move |i| values[*i].clone())
}

pub fn generalized_cigar<T: Hash + Eq + Clone + Display>(
    items: impl Iterator<Item = T>,
    keep_order: bool,
) -> String {
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
                .map(|(item, count)| format!("{}{}", count, item)),
            "",
        )
    }
}

pub fn evidence_kass_raftery_to_letter(evidence: KassRaftery) -> char {
    match evidence {
        KassRaftery::Barely => 'B',
        KassRaftery::None => 'N',
        KassRaftery::Positive => 'P',
        KassRaftery::Strong => 'S',
        KassRaftery::VeryStrong => 'V',
    }
}

/// Collect variants from a given ´bcf::Record`.
pub fn collect_variants(
    record: &mut bcf::Record,
    omit_snvs: bool,
    omit_indels: bool,
    indel_len_range: Option<Range<u64>>,
) -> Result<Vec<model::Variant>> {
    let pos = record.pos() as u64;
    let svlens = match record.info(b"SVLEN").integer() {
        Ok(Some(svlens)) => Some(
            svlens
                .into_iter()
                .map(|l| {
                    if !l.is_missing() {
                        Some(l.abs() as u64)
                    } else {
                        None
                    }
                })
                .collect_vec(),
        ),
        _ => None,
    };
    let end = match record.info(b"END").integer() {
        Ok(Some(end)) => {
            let end = end[0] as u64 - 1;
            Some(end)
        }
        _ => None,
    };
    // TODO avoid cloning svtype
    let svtype = match record.info(b"SVTYPE").string() {
        Ok(Some(svtype)) => Some(svtype[0].to_owned()),
        _ => None,
    };

    // check if len is within the given range
    let is_valid_len = |svlen| {
        if let Some(ref len_range) = indel_len_range {
            // TODO replace with Range::contains once stabilized
            if svlen < len_range.start || svlen >= len_range.end {
                return false;
            }
        }
        true
    };

    let is_valid_insertion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
        alt_allele == b"<INS>"
            || (ref_allele.len() < alt_allele.len()
                && ref_allele == &alt_allele[..ref_allele.len()])
    };

    let is_valid_deletion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
        alt_allele == b"<DEL>"
            || (ref_allele.len() > alt_allele.len()
                && &ref_allele[..alt_allele.len()] == alt_allele)
    };

    let mut variants = Vec::new();

    if let Some(svtype) = svtype {
        if !omit_indels {
            if svtype == b"INS" {
                // get sequence
                let alleles = record.alleles();
                if alleles.len() > 2 {
                    return Err(errors::Error::InvalidBCFRecord {
                        msg: "SVTYPE=INS but more than one ALT allele".to_owned(),
                    })?;
                }
                let ref_allele = alleles[0];
                let alt_allele = alleles[1];

                if alt_allele != b"<INS>" {
                    // don't support insertions without exact sequence
                    let len = alt_allele.len() - ref_allele.len();

                    if is_valid_insertion_alleles(ref_allele, alt_allele)
                        && is_valid_len(len as u64)
                    {
                        variants.push(model::Variant::Insertion(
                            alt_allele[ref_allele.len()..].to_owned(),
                        ));
                    }
                }
            } else if svtype == b"DEL" {
                let svlen = match (svlens, end) {
                    (Some(ref svlens), _) if svlens[0].is_some() => svlens[0].unwrap(),
                    (None, Some(end)) => end - (pos + 1), // pos is pointing to the allele before the DEL
                    _ => {
                        return Err(errors::Error::MissingBCFTag {
                            name: "SVLEN or END".to_owned(),
                        })?;
                    }
                };
                if svlen == 0 {
                    return Err(errors::Error::InvalidBCFRecord {
                        msg: "Absolute value of SVLEN or END - POS must be greater than zero."
                            .to_owned(),
                    })?;
                }
                let alleles = record.alleles();
                if alleles.len() > 2 {
                    return Err(errors::Error::InvalidBCFRecord {
                        msg: "SVTYPE=DEL but more than one ALT allele".to_owned(),
                    })?;
                }
                let ref_allele = alleles[0];
                let alt_allele = alleles[1];

                if alt_allele == b"<DEL>" || is_valid_deletion_alleles(ref_allele, alt_allele) {
                    if is_valid_len(svlen) {
                        variants.push(model::Variant::Deletion(svlen));
                    }
                }
            }
        }
    } else {
        let alleles = record.alleles();
        let ref_allele = alleles[0];

        for (i, alt_allele) in alleles.iter().skip(1).enumerate() {
            if alt_allele == b"<*>" {
                // dummy non-ref allele, signifying potential homozygous reference site
                if !omit_snvs {
                    variants.push(model::Variant::None)
                }
            } else if alt_allele == b"<DEL>" {
                if let Some(ref svlens) = svlens {
                    if let Some(svlen) = svlens[i] {
                        variants.push(model::Variant::Deletion(svlen));
                    }
                    // TODO fail with an error in else case
                }
            } else if alt_allele[0] == b'<' {
                // skip any other special alleles
                ()
            } else if alt_allele.len() == 1 && ref_allele.len() == 1 {
                // SNV
                if !omit_snvs {
                    variants.push(model::Variant::SNV(alt_allele[0]));
                }
            } else if alt_allele.len() == ref_allele.len() {
                // MNV
                variants.push(model::Variant::MNV(alt_allele.to_vec()));
            } else {
                let indel_len = (alt_allele.len() as i64 - ref_allele.len() as i64).abs() as u64;
                // TODO fix position if variant is like this: cttt -> ct

                if omit_indels || !is_valid_len(indel_len) {
                    // skip
                    ()
                } else if is_valid_deletion_alleles(ref_allele, alt_allele) {
                    variants.push(model::Variant::Deletion(
                        (ref_allele.len() - alt_allele.len()) as u64,
                    ));
                } else if is_valid_insertion_alleles(ref_allele, alt_allele) {
                    variants.push(model::Variant::Insertion(
                        alt_allele[ref_allele.len()..].to_owned(),
                    ));
                }
            }
        }
    }

    Ok(variants)
}

/// Sum up in log space the probabilities of the given tags for all variants of
/// vartype in the given BCF record.
///
/// # Arguments
///
/// * `record` - BCF record
/// * `tags` - tags of the set of events to sum up for a particular site and variant
/// * `vartype` - the variant type to consider
pub fn tags_prob_sum(
    record: &mut bcf::Record,
    tags: &[String],
    vartype: Option<&model::VariantType>,
) -> Result<Vec<Option<LogProb>>> {
    let variants = (utils::collect_variants(record, false, false, None))?;
    let mut tags_probs_out = vec![Vec::new(); variants.len()];

    for tag in tags {
        if let Some(tags_probs_in) = (record.info(tag.as_bytes()).float())? {
            //tag present
            for (i, (variant, tag_prob)) in
                variants.iter().zip(tags_probs_in.into_iter()).enumerate()
            {
                if (vartype.is_some() && !variant.is_type(vartype.unwrap())) || tag_prob.is_nan() {
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

pub fn events_to_tags<E>(events: &[E]) -> Vec<String>
where
    E: Event,
{
    events.iter().map(|e| e.tag_name("PROB")).collect_vec()
}

/// Collect distribution of posterior probabilities from a VCF file that has been written by
/// varlociraptor.
///
/// # Arguments
///
/// * `calls` - BCF reader with varlociraptor calls
/// * `events` - the set of events to sum up for a particular site
/// * `vartype` - the variant type to consider
pub fn collect_prob_dist<E>(
    calls: &mut bcf::Reader,
    events: &[E],
    vartype: &model::VariantType,
) -> Result<Vec<NotNan<f64>>>
where
    E: Event,
{
    let mut record = calls.empty_record();
    let mut prob_dist = Vec::new();
    let tags = events_to_tags(events);
    loop {
        if !calls.read(&mut record)? {
            break;
        }

        for p in utils::tags_prob_sum(&mut record, &tags, Some(&vartype))? {
            if let Some(p) = p {
                prob_dist.push(NotNan::new(*p)?);
            }
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
/// * `vartype` - the variant type to consider
pub fn filter_by_threshold<E: Event>(
    calls: &mut bcf::Reader,
    threshold: Option<LogProb>,
    out: &mut bcf::Writer,
    events: &[E],
    vartype: &model::VariantType,
) -> Result<()> {
    let tags = events.iter().map(|e| e.tag_name("PROB")).collect_vec();
    let filter = |record: &mut bcf::Record| {
        let probs = utils::tags_prob_sum(record, &tags, Some(vartype))?;
        Ok(probs.into_iter().map(|p| {
            match (p, threshold) {
                // we allow some numerical instability in case of equality
                (Some(p), Some(threshold)) if p > threshold || relative_eq!(*p, *threshold) => true,
                (Some(_), None) => true,
                _ => false,
            }
        }))
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
pub fn filter_calls<F, I, II>(
    calls: &mut bcf::Reader,
    out: &mut bcf::Writer,
    filter: F,
) -> Result<()>
where
    F: Fn(&mut bcf::Record) -> Result<II>,
    I: Iterator<Item = bool>,
    II: IntoIterator<Item = bool, IntoIter = I>,
{
    let mut record = calls.empty_record();
    loop {
        if !calls.read(&mut record)? {
            return Ok(());
        }

        let mut remove = vec![false]; // don't remove the reference allele
        remove.extend(filter(&mut record)?.into_iter().map(|keep| !keep));

        assert_eq!(
            remove.len(),
            record.allele_count() as usize,
            "bug: filter passed to filter_calls has to return a bool for each alt allele."
        );

        // Write trimmed record if any allele remains. Otherwise skip the record.
        if !remove[1..].iter().all(|r| *r) {
            record.remove_alleles(&remove)?;
            out.write(&record)?;
        }
    }
}

/// Return the greater of two given probabilities.
pub fn max_prob(prob_a: LogProb, prob_b: LogProb) -> LogProb {
    LogProb(*cmp::max(NotNan::from(prob_a), NotNan::from(prob_b)))
}

/// Returns true if all PROB_{event}s are PHRED scaled
pub fn is_phred_scaled(inbcf: &bcf::Reader) -> bool {
    get_event_tags(inbcf)
        .iter()
        // check for missing closing parenthesis for backward compatibility
        .all(|(_, d)| d.ends_with("(PHRED)") || !d.ends_with(")"))
}

/// Returns (ID, Description) for each PROB_{event} INFO tag
pub fn get_event_tags(inbcf: &bcf::Reader) -> Vec<(String, String)> {
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

/// Returns true if given variant is located in a repeat region.
pub fn is_repeat_variant(start: u64, variant: &model::Variant, chrom_seq: &[u8]) -> bool {
    let end = variant.end(start) as usize;
    for nuc in &chrom_seq[start as usize..end] {
        if (*nuc as char).is_lowercase() {
            return true;
        }
    }

    false
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum MiniLogProb {
    F16(f16),
    F32(f32),
}

impl MiniLogProb {
    /// Convert LogProb into a minimal representation for storage.
    /// If integer part is less than -1 and can be represented in f16,
    /// we use f16. Else, we use f32.
    pub fn new(prob: LogProb) -> Self {
        let half = f16::from_f64(*prob);
        let proj = half.to_f64();
        if *prob < -10.0 && proj.floor() as i64 == prob.floor() as i64 {
            MiniLogProb::F16(half)
        } else {
            MiniLogProb::F32(*prob as f32)
        }
    }

    pub fn to_logprob(&self) -> LogProb {
        LogProb(match self {
            MiniLogProb::F16(p) => p.to_f64(),
            MiniLogProb::F32(p) => *p as f64,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::model::VariantType;
    use crate::ComplementEvent;
    use crate::SimpleEvent;
    use bio::stats::{LogProb, Prob};
    use rust_htslib::bcf::{self, Read};

    #[test]
    fn test_tags_prob_sum() {
        // set up test input
        let test_file = "tests/resources/test_tags_prob_sum/overshoot.vcf";
        let mut overshoot_calls = bcf::Reader::from_path(test_file).unwrap();
        let mut record = overshoot_calls.empty_record();
        if let Err(e) = overshoot_calls.read(&mut record) {
            panic!("BCF reading error: {}", e);
        }

        // set up all alt events with names as in prosolo
        let alt_tags = [
            String::from("PROB_ADO_TO_REF"),
            String::from("PROB_ADO_TO_ALT"),
            String::from("PROB_HOM_ALT"),
            String::from("PROB_HET"),
            String::from("PROB_ERR_REF"),
        ];

        let snv = VariantType::SNV;

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
        let absent_event = vec![ComplementEvent {
            name: "absent".to_owned(),
        }];

        let test_file = "tests/resources/test_collect_prob_dist/min.calls.vcf";

        //TESTS deletion
        let del = VariantType::Deletion(None);

        let mut del_calls_1 = bcf::Reader::from_path(test_file).unwrap();
        let prob_del = collect_prob_dist(&mut del_calls_1, &events, &del).unwrap();
        println!("prob_del[0]: {:?}", prob_del[0].into_inner());
        assert_eq!(prob_del.len(), 1);
        assert_relative_eq!(prob_del[0].into_inner(), Prob(0.8).ln(), epsilon = 0.000005);

        let mut del_calls_2 = bcf::Reader::from_path(test_file).unwrap();
        let prob_del_abs = collect_prob_dist(&mut del_calls_2, &absent_event, &del).unwrap();
        assert_eq!(prob_del_abs.len(), 1);
        assert_relative_eq!(
            prob_del_abs[0].into_inner(),
            Prob(0.2).ln(),
            epsilon = 0.000005
        );

        //TESTS insertion
        let ins = VariantType::Insertion(None);

        let mut ins_calls_1 = bcf::Reader::from_path(test_file).unwrap();
        let prob_ins = collect_prob_dist(&mut ins_calls_1, &events, &ins).unwrap();
        assert_eq!(prob_ins.len(), 1);
        assert_relative_eq!(prob_ins[0].into_inner(), Prob(0.2).ln(), epsilon = 0.000005);

        let mut ins_calls_2 = bcf::Reader::from_path(test_file).unwrap();
        let prob_ins_abs = collect_prob_dist(&mut ins_calls_2, &absent_event, &ins).unwrap();
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

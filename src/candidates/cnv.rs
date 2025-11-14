use anyhow::{Context, Result};
use bio_types::genome::{AbstractInterval, AbstractLocus, Interval, Locus};
use itertools::Itertools;
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::{record::Numeric, Format, Header, Read, Reader, Record, Writer};
use rust_htslib::bgzf;
use std::collections::{BTreeMap, HashSet};
use std::ops::Range;
use std::path::PathBuf;
use yaml_rust::yaml::Hash;

use crate::utils::aux_info::AuxInfoCollector;
use crate::variants::model::VariantPrecision;
use crate::variants::types::breakends::Breakend;

#[derive(Builder, Getters, Debug)]
// Computed interval out of the breakends
pub struct BreakendInterval {
    #[getset(get = "pub")]
    interval: Interval,
    #[getset(get = "pub")]
    info_probs: BTreeMap<String, f64>,
    #[getset(get = "pub")]
    involved_afs: Vec<f64>,
}

impl PartialEq for BreakendInterval {
    fn eq(&self, other: &Self) -> bool {
        self.interval == other.interval
            && self.info_probs == other.info_probs
            && self.involved_afs == other.involved_afs
    }
}
impl Eq for BreakendInterval {}

impl std::hash::Hash for BreakendInterval {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.interval.hash(state);
        for v in &self.involved_afs {
            v.to_bits().hash(state); // f64 -> u64 für konsistentes Hashing
        }

        // Hash für info_probs
        for (k, v) in &self.info_probs {
            k.hash(state);
            v.to_bits().hash(state); // f64 -> u64 für konsistentes Hashing
        }
    }
}

#[derive(Builder, Getters, Debug, Clone)]
// Representative for every breakend in the call file
pub struct CnvCandidate {
    #[getset(get = "pub")]
    breakend: Breakend,
    // Map with all relevant info fields as keys and their values in PHRED scale
    #[getset(get = "pub")]
    info_probs: BTreeMap<String, f64>,
    #[getset(get = "pub")]
    af: f64,
}

#[derive(Clone, Debug, Default)]
// For every border, keep track of the accumulated probabilities of all intervals starting/ending there
struct BorderProbs {
    prob_interval_start: f64,
    prob_interval_end: f64,
    // Number of intervals with probability zero
    count_prob_zero: i32,
}

impl BorderProbs {
    fn update_start(&mut self, value: f64) {
        if value != f64::NEG_INFINITY {
            self.prob_interval_start += value;
        } else {
            self.count_prob_zero += 1;
        }
    }

    fn update_end(&mut self, value: f64) {
        if value != f64::NEG_INFINITY {
            self.prob_interval_end += value;
        } else {
            self.count_prob_zero += 1;
        }
    }
}

#[derive(Clone, Debug, Default)]
struct BorderInfo {
    // How many times does an interval start or end here? 1 for start, -1 for end
    delta: i32,
    info_probs_borders: Option<BTreeMap<String, BorderProbs>>,
    afs_start: Option<Vec<f64>>,
    afs_end: Option<Vec<f64>>,
    breakends: Vec<CnvCandidate>,
}

fn write_cnv_records(
    cnv_candidates: HashSet<BreakendInterval>,
    header: &Header,
    outbcf: Option<PathBuf>,
) -> Result<()> {
    let mut bcf_writer = match outbcf {
        Some(path) => Writer::from_path(path, header, true, Format::Bcf)
            .with_context(|| "Error opening BCF writer".to_string())?,
        None => Writer::from_stdout(header, true, Format::Bcf)
            .with_context(|| "Error opening BCF writer".to_string())?,
    };
    for mut candidate in cnv_candidates {
        let mut cnv_record = bcf_writer.empty_record();
        cnv_record.set_pos(candidate.interval.range().start as i64);
        cnv_record.set_qual(f32::missing());
        // TODO: What is the REF allele?
        cnv_record.set_alleles(&[b"N", b"<CNV>"])?;
        cnv_record
            .push_info_integer(b"END", &[candidate.interval.range().end as i32])
            .with_context(|| "Failed to push END info string")?;
        cnv_record
            .push_info_string(b"SVTYPE", &[b"CNV"])
            .with_context(|| "Failed to push SVTYPE info string")?;
        for (k, v) in &candidate.info_probs {
            cnv_record
                .push_info_float(k.as_bytes(), &[*v as f32])
                .with_context(|| format!("Failed to push {} info string", k))?;
        }
        let af_string = if candidate.involved_afs().is_empty() {
            ".".to_string()
        } else {
            candidate
                .involved_afs()
                .iter()
                .map(|x| format!("{:.6}", x)) // oder x.to_string(), falls keine Nachkommastellen
                .collect::<Vec<_>>()
                .join(",")
        };
        let af_bytes = af_string.as_bytes();
        cnv_record.push_format_string(b"AF", &[af_bytes])?;
        bcf_writer
            .write(&cnv_record)
            .with_context(|| "Failed to write BCF record".to_string())?;
    }
    Ok(())
}

fn create_header_from_existing(old_header: &HeaderView) -> anyhow::Result<Header> {
    // kompletten Header aus der View klonen
    let mut header = Header::from_template(old_header);
    header.push_record(
        b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"Ending position of breakend\">",
    );
    header.remove_format(b"AF");
    header.push_record(
    b"##FORMAT=<ID=AF,Number=.,Type=String,Description=\"Allele Frequencies involved in this CNV\">",
    );
    Ok(header)
}

fn parse_breakend(
    record: &Record,
    header: &HeaderView,
    aux_info_collector: &AuxInfoCollector,
) -> Result<Option<Breakend>> {
    let contig = String::from_utf8(
        header
            .rid2name(record.rid().context("Missing rid")?)
            .context("Invalid RID name")?
            .to_vec(),
    )
    .context("Could not create String from rid")?;

    let mut mateid = None;
    if let Ok(Some(values)) = record.info(b"MATEID").string() {
        mateid = values.first().map(|v| v.to_vec());
    }
    let breakend = Breakend::new(
        Locus::new(contig, record.pos() as u64),
        record.alleles()[0],
        record.alleles()[1],
        &record.id(),
        mateid,
        VariantPrecision::Precise,
        aux_info_collector
            .collect(record)
            .context("Auxiliary Info could not be collected from record")?,
    )
    .context("Breakend could not be created")?;
    Ok(breakend)
}

fn breakends_to_interval_borders(
    candidates: Vec<CnvCandidate>,
    number_combo_intervals: usize,
    // In the first step we do not want to apply all probs to the borders since we do not autoinclude all candidates over this border.
    // In the second step we do want to apply all probs to the borders since we want to compute the final intervals for specific combinations of candidates
    apply_probs: bool,
) -> BTreeMap<Locus, BorderInfo> {
    let mut all_intervals: BTreeMap<Locus, BorderInfo> = BTreeMap::new();
    for c in candidates {
        // Find left and right side of the interval
        if let Some(join) = c.breakend.join() {
            let (bp_start, bp_end) = (c.breakend.locus(), join.locus());
            let (interval_start, interval_end) = if bp_start.pos() < bp_end.pos() {
                (bp_start.to_owned(), bp_end.to_owned())
            } else {
                (bp_end.to_owned(), bp_start.to_owned())
            };

            // Start-Border, whenever two intervals start or end at the same position, accumulate their probabilities
            let entry_start = all_intervals.entry(interval_start).or_default();
            entry_start.delta += 1;
            // if number_combo_intervals == 1 {
            match &mut entry_start.afs_start {
                Some(afs) => afs.push(c.af),
                None => entry_start.afs_start = Some(vec![c.af]),
            }
            // }

            if apply_probs {
                let map = entry_start
                    .info_probs_borders
                    .get_or_insert_with(BTreeMap::new);

                for (k, &v) in &c.info_probs {
                    map.entry(k.clone()).or_default().update_start(v);
                }
            }
            entry_start.breakends.push(c.clone());

            // End-Border, whenever two intervals start or end at the same position, accumulate their probabilities
            let entry_end = all_intervals.entry(interval_end).or_default();
            entry_end.delta -= 1;
            // if number_combo_intervals == 1 {
            match &mut entry_end.afs_end {
                Some(afs) => afs.push(c.af),
                None => entry_end.afs_end = Some(vec![c.af]),
            }
            // }

            // entry_end.breakends.push(&c);
            if apply_probs {
                let map = entry_end
                    .info_probs_borders
                    .get_or_insert_with(BTreeMap::new);

                for (k, v) in c.clone().info_probs {
                    map.entry(k.clone()).or_default().update_end(v);
                }
            }
        }
    }
    all_intervals
}

fn compute_overlapping_intervals(candidates: &Vec<CnvCandidate>) -> HashSet<BreakendInterval> {
    let mut common_intervals = HashSet::new();
    for k in 1..=candidates.len() {
        for combo in candidates.iter().combinations(k) {
            let combo_owned: Vec<CnvCandidate> = combo.into_iter().cloned().collect();
            let interval_borders = breakends_to_interval_borders(combo_owned, k, true);
            let new_intervals = interval_borders_to_intervals_with_probs(interval_borders);
            common_intervals.extend(new_intervals);
        }
    }
    common_intervals
}

fn interval_borders_to_intervals_with_probs(
    borders: BTreeMap<Locus, BorderInfo>,
) -> Vec<BreakendInterval> {
    let mut cnv_intervals = Vec::new();
    let mut last_pos: Option<Locus> = None;
    let mut interval_cov = 0;
    let mut interval_info_probs: BTreeMap<String, f64> = BTreeMap::new();
    let mut border_afs: Vec<f64> = Vec::new();
    for (locus, border_info) in &borders {
        // If we have a previous position and the coverage is not zero, create an interval
        if let Some(prev) = &last_pos {
            if prev.contig() == locus.contig() && prev.pos() < locus.pos() && interval_cov != 0 {
                cnv_intervals.push(BreakendInterval {
                    interval: Interval::new(
                        prev.contig().to_owned(),
                        Range {
                            start: prev.pos() + 1,
                            end: locus.pos(),
                        },
                    ),
                    info_probs: interval_info_probs.clone(),
                    involved_afs: border_afs.clone(),
                });
                dbg!(&cnv_intervals.last());
            }
        }

        // Update current_info_probs by adding the start probabilities and
        // subtracting the end probabilities at this border
        for (k, v) in &border_info.info_probs_borders.clone().unwrap() {
            interval_info_probs
                .entry(k.clone())
                .and_modify(|prob| {
                    // If any of the intervals contributing to this border has probability zero, the resulting interval has probability zero
                    *prob = if v.count_prob_zero > 0 {
                        0.0
                    } else {
                        let delta = v.prob_interval_start - v.prob_interval_end;
                        let delta = if delta.is_finite() { delta } else { 0.0 };

                        *prob + delta
                    };
                })
                .or_insert_with(|| {
                    if v.count_prob_zero > 0 {
                        0.0
                    } else {
                        v.prob_interval_start - v.prob_interval_end
                    }
                });
        }
        border_afs.extend(border_info.afs_start.clone().unwrap_or_else(|| vec![]));

        dbg!(&border_afs);
        let to_remove = border_info.afs_end.clone().unwrap_or_else(Vec::new);
        border_afs.retain(|x| !to_remove.contains(x));
        interval_cov += border_info.delta;
        last_pos = Some(locus.clone());
    }
    cnv_intervals
}

fn borders_to_intervals(
    borders: BTreeMap<Locus, BorderInfo>,
    apply_probs: bool,
) -> HashSet<BreakendInterval> {
    let mut cnv_intervals: HashSet<BreakendInterval> = HashSet::new();
    let mut current_candidates = Vec::new();
    let mut last_pos: Option<Locus> = None;
    let mut interval_cov = 0;
    let mut interval_info_probs: BTreeMap<String, f64> = BTreeMap::new();

    for (locus, border_info) in &borders {
        // If we have a previous position and the coverage is not zero, create an interval
        interval_cov += border_info.delta;
        if let Some(prev) = &last_pos {
            if prev.contig() == locus.contig() && prev.pos() < locus.pos() {
                // Apply probability adjustments
                if interval_cov == 0 {
                    // cnv_intervals.push(BreakendInterval {
                    //     interval: Interval::new(
                    //         prev.contig().to_owned(),
                    //         Range {
                    //             start: prev.pos() + 1,
                    //             end: locus.pos(),
                    //         },
                    //     ),
                    //     info_probs: interval_info_probs.clone(),
                    // });
                    cnv_intervals.extend(compute_overlapping_intervals(&current_candidates));
                    current_candidates.clear();

                    // TODO: Compute common intervals
                } else {
                    current_candidates.extend(border_info.breakends.clone());
                }
            } else {
                warn!("There is a change of contigs or the positions are sorted badly");
            }
        } else {
            current_candidates.extend(border_info.breakends.clone());
        }

        // // Update current_info_probs by adding the start probabilities and
        // // subtracting the end probabilities at this border
        // for (k, v) in &border_info.info_probs_borders {
        //     interval_info_probs
        //         .entry(k.clone())
        //         .and_modify(|prob| {
        //             // If any of the intervals contributing to this border has probability zero, the resulting interval has probability zero
        //             *prob = if v.count_prob_zero > 0 {
        //                 0.0
        //             } else {
        //                 let delta = v.prob_interval_start - v.prob_interval_end;
        //                 let delta = if delta.is_finite() { delta } else { 0.0 };

        //                 *prob + delta
        //             };
        //         })
        //         .or_insert_with(|| {
        //             if v.count_prob_zero > 0 {
        //                 0.0
        //             } else {
        //                 v.prob_interval_start - v.prob_interval_end
        //             }
        //         });
        // }
        last_pos = Some(locus.clone());
    }

    cnv_intervals
}

/// Reads an input BCF file, processes each record, and writes out modified CNV records
pub fn find_candidates(breakends_bcf: PathBuf, outbcf: Option<PathBuf>) -> Result<()> {
    let mut bcf_reader = Reader::from_path(&breakends_bcf)
        .with_context(|| format!("Error opening input file: {:?}", breakends_bcf))?;

    // Create header for interval bcf
    let header = create_header_from_existing(bcf_reader.header())?;
    let header_orig = bcf_reader.header().clone();

    let mut breakends = Vec::new();

    // TODO: Compute dynamicall (problem: Cant extract INFO information from Header(View)
    let prob_keys = vec![
        "PROB_HIGH".to_string(),
        "PROB_LOW".to_string(),
        "PROB_ABSENT".to_string(),
        "PROB_ARTIFACT".to_string(),
    ];
    let aux_info = AuxInfoCollector::new(&[], &bcf_reader)
        .with_context(|| "Failed to create AuxInfoCollector")?;

    // Exract breakends + information out of every record in this for loop
    for record_result in bcf_reader.records() {
        let record = record_result.with_context(|| "Error reading record")?;

        let breakend = match parse_breakend(&record, &header_orig, &aux_info)? {
            Some(b) if b.is_left_to_right() && b.join().is_some() => b,
            _ => continue,
        };

        let af = {
            let af_values = record.format(b"AF").float()?;
            if af_values.is_empty() || af_values[0].is_empty() {
                f64::NAN
            } else {
                af_values[0][0] as f64
            }
        };

        let mut probs = BTreeMap::new();
        for key in &prob_keys {
            if let Ok(Some(value)) = record.info(key.as_bytes()).float() {
                if let Some(&val) = value.first() {
                    // Probs in the PROB tags are PHRED encoded, convert to linear scale

                    probs.insert(key.clone(), val as f64);
                }
            }
        }
        // Create new breakend candidate out of record
        let candidate = CnvCandidateBuilder::default()
            .breakend(breakend)
            .info_probs(probs)
            .af(af)
            .build()?;

        breakends.push(candidate);
    }
    let interval_borders = breakends_to_interval_borders(breakends, 0, false);
    // Compute intervals out of all breakends
    let intervals: HashSet<BreakendInterval> = borders_to_intervals(interval_borders, false);
    write_cnv_records(intervals, &header, outbcf)?;
    Ok(())
}

use std::error::Error;
use std::fs;
use std::path::Path;

use itertools::zip;
use itertools::Itertools;
use itertools_num::ItertoolsNum;
use rust_htslib::bcf;
use rust_htslib::bcf::header::HeaderRecord;
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::{Read, Record};
use serde::{Deserialize, Serialize};
use serde_json::json;
use serde_json::value::Value::Object;
use std::collections::HashSet;

#[derive(Serialize, Deserialize)]
struct CNRecord {
    rid: u32,
    start: u32,
    end: u32,
    copynumber: u32,
    subclone_fraction: f32,
}

#[derive(Serialize, Deserialize)]
struct Contig {
    id: String,
    length: usize,
    cumulative_position: usize,
}

/// Plot copynumber vs genomic position (in bp), with subclone fraction as color dimension.
///
/// Uses the vega-lite v3 json scheme.
/// Currently, chromosomes are divided using vertical lines, with genomic positions w.r.t hg38.
pub fn plot<P: AsRef<Path>>(cnv_calls: Option<P>) -> Result<(), Box<Error>> {
    let mut inbcf_reader = match cnv_calls {
        Some(p) => bcf::Reader::from_path(p)?,
        None => bcf::Reader::from_stdin()?,
    };

    let mut cumulative_pos = 0;
    let mut last_record: Option<Record> = None;
    let mut records: Vec<CNRecord> = Vec::new();
    let header = inbcf_reader.header();

    // TODO use user defined ids (for now default to chr1 - chr22, i.e. skip X, Y, M and misc. scaffolds)
    let ids: HashSet<usize> = (0..22).collect();
    let mut contigs = get_contig_info(header, &ids);

    for record in inbcf_reader.records() {
        let mut record = record?;
        if let Some(mut previous_record) = last_record {
            if record.rid() != previous_record.rid() {
                cumulative_pos += previous_record.pos();
            } else {
                records.push(CNRecord {
                    // TODO handle this and other expects below as errors once we have moved to snafu.
                    rid: previous_record.rid().expect("Failed reading rid"),
                    start: previous_record
                        .info(b"END")
                        .integer()?
                        .expect("Failed reading info tag 'END'")[0]
                        as u32
                        + cumulative_pos,
                    end: record.pos() + cumulative_pos,
                    copynumber: 2,
                    subclone_fraction: 0f32,
                });
            }
        }
        let copynumber = record
            .info(b"CN")
            .integer()?
            .expect("Failed reading info tag 'CN'")[0];
        let subclone_fraction = record
            .info(b"VAF")
            .float()?
            .expect("Failed reading info tag 'VAF'")[0];
        let start = record.pos();
        let end = record
            .info(b"END")
            .integer()?
            .expect("Failed reading info tag 'END'")[0] as u32;
        records.push(CNRecord {
            rid: record.rid().expect("Failed reading rid"),
            start: start + cumulative_pos,
            end: end + cumulative_pos,
            copynumber: copynumber as u32,
            subclone_fraction,
        });
        last_record = Some(record);
    }

    let mut cumulative_chr_pos = contigs
        .iter()
        .map(|(_, length, _)| length)
        .cumsum::<usize>()
        .collect_vec();
    cumulative_chr_pos.insert(0, 0);
    contigs.push(("".to_string(), 0, 0));
    let contigs_info = zip(contigs, cumulative_chr_pos)
        .map(|((id, length, _), pos)| Contig {
            id,
            length,
            cumulative_position: pos,
        })
        .collect_vec();

    let records = records.iter().filter(|&r| r.rid < 22).collect_vec();

    // read vega-lite plot json and replace data entry for "copynumbers" with actual values.
    let mut blueprint =
        serde_json::from_str(&fs::read_to_string(&"templates/blueprint_plot_cnv.json")?)?;
    if let Object(ref mut blueprint) = blueprint {
        let datasets = &mut blueprint["datasets"];
        if let Object(ref mut datasets) = datasets {
            datasets.insert("copynumbers".to_owned(), json!(records));
            datasets.insert("contigs".to_owned(), json!(contigs_info));
        }
        println!("{}", serde_json::to_string_pretty(blueprint)?);
    }
    Ok(())
}

fn get_contig_info(header: &HeaderView, ids: &HashSet<usize>) -> Vec<(String, usize, usize)> {
    header
        .header_records()
        .iter()
        .filter_map(|hr| match hr {
            HeaderRecord::Contig { key, values } => {
                let idx = values["IDX"].parse::<usize>().expect("Failed parsing IDX"); // TODO snafu
                if !ids.contains(&idx) {
                    None
                } else {
                    Some((
                        values["ID"].clone(),
                        values["length"]
                            .parse::<usize>()
                            .expect("Failed parsing length"), // TODO snafu
                        idx,
                    ))
                }
            }
            _ => None,
        })
        .collect_vec()
}

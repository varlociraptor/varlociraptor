use std::error::Error;
use std::fs;
use std::path::Path;

use rust_htslib::bcf;
use rust_htslib::bcf::{Read, Record};
use serde::{Deserialize, Serialize};
use serde_json::json;
use serde_json::value::Value::Object;

#[derive(Serialize, Deserialize)]
struct CNRecord {
    start: u32,
    end: u32,
    copynumber: u32,
    subclone_fraction: f32,
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

    for record in inbcf_reader.records() {
        let mut record = record?;
        if let Some(mut previous_record) = last_record {
            if record.rid() != previous_record.rid() {
                cumulative_pos += previous_record.pos();
            } else {
                records.push(CNRecord {
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
            start: start + cumulative_pos,
            end: end + cumulative_pos,
            copynumber: copynumber as u32,
            subclone_fraction,
        });
        last_record = Some(record);
    }

    // read vega-lite plot json and replace data entry for "copynumbers" with actual values.
    let mut blueprint =
        serde_json::from_str(&fs::read_to_string(&"templates/blueprint_plot_cnv.json")?)?;
    if let Object(ref mut blueprint) = blueprint {
        let datasets = &mut blueprint["datasets"];
        if let Object(ref mut datasets) = datasets {
            datasets.insert("copynumbers".to_owned(), json!(records));
        }
        println!("{}", serde_json::to_string_pretty(blueprint)?);
    }
    Ok(())
}

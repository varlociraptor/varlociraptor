use std::error::Error;

use bio::stats::{PHREDProb, Prob};
use itertools::Itertools;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;

/// Decode PHRED scaled values to probabilities.
pub fn decode_phred() -> Result<(), Box<Error>> {
    let mut inbcf = bcf::Reader::from_stdin()?;
    let tags = inbcf
        .header()
        .header_records()
        .into_iter()
        .filter_map(|rec| match rec {
            bcf::header::HeaderRecord::Info { ref key, .. } if key.starts_with("PROB_") => {
                Some(key.clone().into_bytes())
            }
            _ => None,
        })
        .collect_vec();
    // setup output file
    let header = bcf::Header::from_template(inbcf.header());
    let mut outbcf = bcf::Writer::from_stdout(&header, false, false)?;

    for record in inbcf.records() {
        let mut record = record?;
        for tag in &tags {
            if let Some(values) = record.info(tag).float()? {
                let converted = values
                    .iter()
                    .map(|v| *Prob::from(PHREDProb(*v as f64)) as f32)
                    .collect_vec();
                record.clear_info_float(tag)?;
                record.push_info_float(tag, &converted)?;
            }
        }
        outbcf.write(&record)?;
    }

    Ok(())
}

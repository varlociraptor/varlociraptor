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
        .filter_map(|rec| {
            if let bcf::header::HeaderRecord::Info { values, .. } = rec {
                let id = values.get("ID").unwrap().clone();
                if id.starts_with("PROB_") {
                    let description = values.get("Description").unwrap().clone();
                    // only consider PHRED scaled probabilities
                    if description.ends_with("(PHRED)") {
                        return Some((id, description));
                    }
                }
            }
            None
        })
        .collect_vec();

    // setup output file
    let mut header = bcf::Header::from_template(inbcf.header());
    for (id, description) in &tags {
        header.remove_info(&id.clone().into_bytes());
        let description = description.replace("(PHRED)", "(linear)");
        header.push_record(
            format!(
                "##INFO=<ID={},Number=A,Type=Float,\
                 Description=\"{}\">",
                id, description
            )
            .as_bytes(),
        );
    }
    let mut outbcf = bcf::Writer::from_stdout(&header, false, false)?;

    for record in inbcf.records() {
        let mut record = record?;
        for (tag, _desc) in &tags {
            let id = &tag.clone().into_bytes();
            if let Some(values) = record.info(id).float()? {
                let converted = values
                    .iter()
                    .map(|v| *Prob::from(PHREDProb(*v as f64)) as f32)
                    .collect_vec();
                record.push_info_float(&id, &converted)?;
            }
        }
        outbcf.write(&record)?;
    }

    Ok(())
}

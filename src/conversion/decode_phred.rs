use std::error::Error;

use crate::utils::get_event_tags;
use bio::stats::{PHREDProb, Prob};
use itertools::Itertools;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;

/// Decode PHRED scaled values to probabilities.
pub fn decode_phred() -> Result<(), Box<Error>> {
    let mut inbcf = bcf::Reader::from_stdin()?;

    let tags = get_event_tags(&inbcf)
        .iter()
        .filter(|(tag, desc)| desc.ends_with("(PHRED)"))
        .cloned()
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

use std::path::Path;
use std::error::Error;

use bio::stats::{LogProb, PHREDProb, bayesian};
use itertools::Itertools;
use rust_htslib::bcf;

use Event;
use BCFError;


/// Annotate a given VCF/BCF file with FDR estimates for the given events.
///
/// `inbcf` - path to VCF/BCF file that shall be annotated (`"-"` for STDIN)
/// `outbcf` - path to resulting BCF file (`"-"` for STDOUT)
/// `events` - set of events to annotate
pub fn annotate<R, W, E>(inbcf: &R, outbcf: &W, events: &[E]) -> Result<(), Box<Error>> where
    R: AsRef<Path>,
    W: AsRef<Path>,
    E: Event
{
    // read probs
    let reader = try!(bcf::Reader::new(&inbcf));

    let mut event_peps = Vec::new();
    for _ in 0..events.len() {
        event_peps.push(Vec::new());
    }

    let mut record = bcf::Record::new();

    // collect posterior error probabilities
    loop {
        if let Err(e) = reader.read(&mut record) {
            if e.is_eof() {
                break;
            } else {
                return Err(Box::new(e));
            }
        }

        for (event, peps) in events.iter().zip(event_peps.iter_mut()) {
            let tag_name = event.tag_name("PROB");
            let rec_peps = try!(record.info(tag_name.as_bytes()).float());
            match rec_peps {
                Some(rec_peps) => {
                    for &p in rec_peps.iter() {
                        // posterior error probability
                        let pep = LogProb::from(PHREDProb(p as f64)).ln_one_minus_exp();
                        peps.push(pep);
                    }
                },
                None => {
                    return Err(Box::new(BCFError::MissingTag(tag_name.clone())));
                }
            }

        }
    }

    // estimate FDR
    let event_fdrs = event_peps.iter().map(|peps| {
        bayesian::expected_fdr(peps).iter().map(|&p| *PHREDProb::from(p) as f32).collect_vec()
    }).collect_vec();

    // write results

    // read bcf again
    let reader = try!(bcf::Reader::new(&inbcf));
    let mut outbcf = {
        let mut header = bcf::Header::with_template(&reader.header);
        for event in events {
            header.push_record(
                event.header_entry("FDR", "PHRED-scaled FDR when considering this and all better").as_bytes()
            );
        }
        try!(bcf::Writer::new(outbcf, &header, false, false))
    };

    let mut i = 0;
    loop {
        if let Err(e) = reader.read(&mut record) {
            if e.is_eof() {
                break;
            } else {
                return Err(Box::new(e));
            }
        }
        // translate to new header
        outbcf.translate(&mut record);

        let allele_count = record.alleles().len() - 1;

        // write FDRs for all events
        for (event, fdrs) in events.iter().zip(event_fdrs.iter()) {
            try!(record.push_info_float(event.tag_name("FDR").as_bytes(), &fdrs[i..i + allele_count]));
        }

        i += allele_count;
    }

    Ok(())
}

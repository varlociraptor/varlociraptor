use std::ops::Range;
use std::error::Error;
use std::fs;
use std::str;

use itertools::Itertools;
use rust_htslib::bcf;
use bio::io::fasta;
use bio::stats::{LogProb, PHREDProb};
use ordered_float::NotNaN;

use model;
use BCFError;
use utils;
use Event;


/// Collect variants from a given ´bcf::Record`.
pub fn collect_variants(
    record: &mut bcf::Record,
    omit_snvs: bool,
    omit_indels: bool,
    indel_len_range: Option<Range<u32>>,
    exclusive_end: bool
) -> Result<Vec<Option<model::Variant>>, Box<Error>> {
    let pos = record.pos();
    let svlen = match record.info(b"SVLEN").integer() {
        Ok(Some(svlen)) => Some(svlen[0].abs() as u32),
        _ => None
    };
    let end = match record.info(b"END").integer() {
        Ok(Some(end)) => {
            let mut end = end[0] as u32 - 1;
            if exclusive_end {
                // this happens with DELLY
                debug!("fixing END tag");
                end -= 1;
            }
            Some(end)
        },
        _ => None
    };
    // TODO avoid cloning svtype
    let svtype = match record.info(b"SVTYPE").string() {
        Ok(Some(svtype)) => Some(svtype[0].to_owned()),
        _ => None
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

    let variants = if let Some(svtype) = svtype {
        vec![
            if omit_indels {
                None
            } else if svtype == b"INS" {
                // get sequence
                let alleles = record.alleles();
                if alleles.len() > 2 {
                    return Err(Box::new(
                        BCFError::InvalidRecord("SVTYPE=INS but more than one ALT allele".to_owned())
                    ))
                }
                let ref_allele = alleles[0];
                let alt_allele = alleles[1];

                if alt_allele == b"<INS>" {
                    // don't support insertions without exact sequence
                    None
                } else {
                    let len = alt_allele.len() - ref_allele.len();

                    if is_valid_len(len as u32) {
                        Some(model::Variant::Insertion(alt_allele[ref_allele.len()..].to_owned()))
                    } else {
                        None
                    }
                }
            } else if svtype == b"DEL" {
                let svlen = match(svlen, end) {
                    (Some(svlen), _)  => svlen,
                    (None, Some(end)) => end - pos,
                    _ => {
                        return Err(Box::new(BCFError::MissingTag("SVLEN or END".to_owned())));
                    }
                };
                if is_valid_len(svlen) {
                    Some(model::Variant::Deletion(svlen))
                } else {
                    None
                }
            } else {
                None
            }
        ]
    } else {
        let alleles = record.alleles();
        let ref_allele = alleles[0];

        alleles.iter().skip(1).map(|alt_allele| {
            if alt_allele[0] == b'<' {
                // skip allele if it is a special tag (such alleles have been handled above)
                None
            } else if alt_allele.len() == 1 && ref_allele.len() == 1 {
                // SNV
                if omit_snvs {
                    None
                } else {
                    Some(model::Variant::SNV(alt_allele[0]))
                }
            } else if alt_allele.len() == ref_allele.len() {
                // neither indel nor SNV
                None
            } else {
                let indel_len = (alt_allele.len() as i32 - ref_allele.len() as i32).abs() as u32;

                if omit_indels {
                    None
                } else if !is_valid_len(indel_len) {
                    None
                } else if alt_allele.len() < ref_allele.len() {
                    Some(model::Variant::Deletion((ref_allele.len() - alt_allele.len()) as u32))
                } else {
                    Some(model::Variant::Insertion(alt_allele[ref_allele.len()..].to_owned()))
                }
            }
        }).collect_vec()
    };

    Ok(variants)
}


/// A lazy buffer for reference sequences.
pub struct ReferenceBuffer {
    reader: fasta::IndexedReader<fs::File>,
    chrom: Option<Vec<u8>>,
    sequence: Vec<u8>
}


impl ReferenceBuffer {
    pub fn new(fasta: fasta::IndexedReader<fs::File>) -> Self {
        ReferenceBuffer {
            reader: fasta,
            chrom: None,
            sequence: Vec::new()
        }
    }

    /// Load given chromosome and return it as a slice. This is O(1) if chromosome was loaded before.
    pub fn seq(&mut self, chrom: &[u8]) -> Result<&[u8], Box<Error>> {
        if let Some(ref last_chrom) = self.chrom {
            if last_chrom == &chrom {
                return Ok(&self.sequence);
            }
        }

        try!(self.reader.read_all(try!(str::from_utf8(chrom)), &mut self.sequence));
        self.chrom = Some(chrom.to_owned());

        Ok(&self.sequence)
    }
}


/// Collect distribution of posterior probabilities from a VCF file that has been written by
/// libprosic.
///
/// # Arguments
///
/// * `calls` - BCF reader with libprosic calls
/// * `events` - the set of events to sum up for a particular site
/// * `vartype` - the variant type to consider
pub fn collect_prob_dist<E: Event>(
    calls: &mut bcf::Reader,
    events: &[E],
    vartype: &model::VariantType) -> Result<Vec<NotNaN<f64>>, Box<Error>> {
    let mut record = bcf::Record::new();
    let mut prob_dist = Vec::new();
    let tags = events.iter().map(|e| e.tag_name("PROB")).collect_vec();
    loop {
        if let Err(e) = calls.read(&mut record) {
            if e.is_eof() {
                break;
            } else {
                return Err(Box::new(e));
            }
        }

        let variants = (utils::collect_variants(&mut record, false, false, None, false))?;
        let mut events_prob_sum = LogProb::ln_zero();
        for tag in &tags {
            if let Some(event_probs) = (record.info(tag.as_bytes()).float())? {
                //tag present
                for (variant, event_prob) in (&variants).into_iter().zip(event_probs.into_iter()) {
                    if let Some(ref variant) = *variant {
                        if !variant.is_type(vartype) || event_prob.is_nan() {
                            continue;
                        }
                        events_prob_sum = events_prob_sum.ln_add_exp( LogProb::from( PHREDProb( *event_prob as f64 ) ) );
                    }
                }

            }
        }
        prob_dist.push(try!(NotNaN::new( *events_prob_sum ) ));
    }
    prob_dist.sort();
    Ok(prob_dist)
}

/// Filter a VCF record stream by a minimum threshold on the sum of
/// posterior probabilities of a given set of Events. The threshold
/// should be an informative false discovery rate (FDR) threshold,
/// e.g. determined with the libprosic FDR control functionality.
///
/// # Arguments
///
/// * `calls` - BCF reader with libprosic calls
/// * `threshold` - minimum threshold for the sum of posterior probabilities of the set of Events considered
/// * `calls` - BCF writer for the filtered libprosic calls
/// * `events` - the set of Events to filter on
/// * `vartype` - the variant type to consider
pub fn filter_by_threshold<E: Event>(
    calls: &mut bcf::Reader,
    threshold: &f64,
    out: &mut bcf::Writer,
    events: &[E],
    vartype: &model::VariantType
) -> Result<(), Box<Error>> {
    let mut record = bcf::Record::new();
    let tags = events.iter().map(|e| e.tag_name("PROB")).collect_vec();
    let lp_threshold = LogProb::from( PHREDProb( *threshold + 0.000000001) ); // the manual epsilon is required, because the threshold output by `control-fdr` has some digits cut off, which can lead to the threshold being lower than the values reread from the BCF record only due to a higher precision
    loop {
        if let Err(e) = calls.read(&mut record) {
            if e.is_eof() {
                return Ok(());
            } else {
                return Err(Box::new(e));
            }
        }

        let variants = (utils::collect_variants(&mut record, false, false, None, false))?;
        let mut events_probs = Vec::with_capacity( variants.len() * tags.len() );
        for tag in &tags {
            if let Some(event_probs) = (record.info(tag.as_bytes()).float())? {
                //tag present
                for (variant, event_prob) in (&variants).into_iter().zip(event_probs.into_iter()) {
                    if let Some(ref variant) = *variant {
                        if !variant.is_type(vartype) || event_prob.is_nan() {
                            continue;
                        }
                        events_probs.push( LogProb::from( PHREDProb( *event_prob as f64 ) ) );
                    }
                }
            }
        }
        let events_prob_sum = LogProb::ln_sum_exp(&events_probs);
        if events_prob_sum >= lp_threshold { out.write(&record)? };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use rust_htslib::bcf;
    use bio::stats::Prob;
    use model::VariantType;
    use ComplementEvent;
    use SimpleEvent;

    #[test]
    fn test_collect_prob_dist() {
        // setup events with names as in prosic2
        let events = vec![
            SimpleEvent { name: "germline".to_owned() },
            SimpleEvent { name: "somatic".to_owned() }
        ];
        // setup absent event as the complement of the other events
        let absent_event = vec![ ComplementEvent { name: "absent".to_owned() } ];

        let test_file = "tests/resources/test_collect_prob_dist/min.calls.vcf";

        //TESTS deletion
        let del = VariantType::Deletion(None);

        let mut del_calls_1 = bcf::Reader::from_path( test_file ).unwrap();
        if let Ok(prob_del) = collect_prob_dist(&mut del_calls_1, &events, &del) {
            println!("prob_del[0]: {:?}", prob_del[0].into_inner() );
            assert_eq!( prob_del.len(), 3 );
            assert_relative_eq!( prob_del[2].into_inner(), Prob(0.8).ln(), epsilon = 0.000005 );
        } else {
            panic!("collect_prob_dist(&calls, &events, &del) returned Error")
        }
        let mut del_calls_2 = bcf::Reader::from_path( test_file ).unwrap();
        if let Ok(prob_del_abs) = collect_prob_dist(&mut del_calls_2, &absent_event, &del) {
            assert_eq!( prob_del_abs.len(), 3 );
            assert_relative_eq!( prob_del_abs[2].into_inner(), Prob(0.2).ln(), epsilon = 0.000005 );
        } else {
            panic!("collect_prob_dist(&calls, &absent_event, &del) returned Error")
        }

        //TESTS insertion
        let ins = VariantType::Insertion(None);

        let mut ins_calls_1 = bcf::Reader::from_path( test_file ).unwrap();
        if let Ok(prob_ins) = collect_prob_dist(&mut ins_calls_1, &events, &ins) {
            assert_eq!( prob_ins.len(), 3 );
            assert_relative_eq!( prob_ins[2].into_inner(), Prob(0.2).ln(), epsilon = 0.000005 );
        } else {
            panic!("collect_prob_dist(&calls, &events, &ins) returned Error")
        }
        let mut ins_calls_2 = bcf::Reader::from_path( test_file ).unwrap();
        if let Ok(prob_ins_abs) = collect_prob_dist(&mut ins_calls_2, &absent_event, &ins) {
            assert_eq!( prob_ins_abs.len(), 3 );
            assert_relative_eq!( prob_ins_abs[2].into_inner(), Prob(0.8).ln(), epsilon = 0.000005 );
        } else {
            panic!("collect_prob_dist(&calls, &absent_event, &ins) returned Error")
        }
    }

}

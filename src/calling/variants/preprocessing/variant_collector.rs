use bio_types::genome;

use crate::variants::{Variant, self};

pub struct VariantCollector {
    omit_snvs: bool,
    omit_indels: bool,
    reference_buffer: ReferenceBuffer,
    realigner: 
}

impl VariantCollector {
    pub fn collect(&mut self, contig_name: &[u8], ref_seq: &[u8], realigner: ) -> Vec<Box<dyn Variant>> {
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

        let locus = || {
            genome::Locus::new(contig_name.to_owned(), pos)
        };

        let interval = |end: u64| {
            genome::Interval::new(contig_name.to_owned(), pos, end)
        }

        let variants = if let Some(svtype) = svtype {
            vec![if omit_indels {
                None
            } else if svtype == b"INS" {
                // get sequence
                let alleles = record.alleles();
                if alleles.len() > 2 {
                    return Err(errors::Error::InvalidBCFRecord {
                        msg: "SVTYPE=INS but more than one ALT allele".to_owned(),
                    })?;
                }
                let ref_allele = alleles[0];
                let alt_allele = alleles[1];

                if alt_allele == b"<INS>" {
                    // don't support insertions without exact sequence
                    None
                } else {
                    let len = alt_allele.len() - ref_allele.len();

                    if is_valid_insertion_alleles(ref_allele, alt_allele) && is_valid_len(len as u64) {
                        Some(model::Variant::Insertion(
                            alt_allele[ref_allele.len()..].to_owned(),
                        ))
                    } else {
                        None
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
                        Some(model::Variant::Deletion(svlen))
                    } else {
                        None
                    }
                } else {
                    None
                }
            } else {
                None
            }]
        } else {
            let alleles = record.alleles();
            let ref_allele = alleles[0];

            alleles
                .iter()
                .skip(1)
                .enumerate()
                .filter_map(|(i, alt_allele)| {
                    if alt_allele == b"<*>" {
                        // dummy non-ref allele, signifying potential homozygous reference site
                        if omit_snvs {
                            None
                        } else {
                            Some(Box::new(variants::None::new(locus(), ref_seq[start])))
                        }
                    } else if alt_allele == b"<DEL>" {
                        if let Some(ref svlens) = svlens {
                            if let Some(svlen) = svlens[i] {
                                Some(Box::new(variants::Deletion::new(interval(), realigner)))
                            } else {
                                // TODO fail with an error in this case
                                None
                            }
                        } else {
                            // TODO fail with an error in this case
                            None
                        }
                    } else if alt_allele[0] == b'<' {
                        // skip any other special alleles
                        None
                    } else if alt_allele.len() == 1 && ref_allele.len() == 1 {
                        // SNV
                        if omit_snvs {
                            None
                        } else {
                            Some(model::Variant::SNV(alt_allele[0]))
                        }
                    } else if alt_allele.len() == ref_allele.len() {
                        // MNV
                        Some(model::Variant::MNV(alt_allele.to_vec()))
                    } else {
                        let indel_len =
                            (alt_allele.len() as i64 - ref_allele.len() as i64).abs() as u64;
                        // TODO fix position if variant is like this: cttt -> ct

                        if omit_indels || !is_valid_len(indel_len) {
                            None
                        } else if is_valid_deletion_alleles(ref_allele, alt_allele) {
                            Some(model::Variant::Deletion(
                                (ref_allele.len() - alt_allele.len()) as u64,
                            ))
                        } else if is_valid_insertion_alleles(ref_allele, alt_allele) {
                            Some(model::Variant::Insertion(
                                alt_allele[ref_allele.len()..].to_owned(),
                            ))
                        } else {
                            None
                        }
                    }
                })
                .collect_vec()
        };

        Ok(variants)
    }
}
//! bcf_utils.rs
//!
//! Utilities for handling/parsing BCF files analysis.

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, Result};
use log::{info, warn};
use rust_htslib::bcf::{self, header::HeaderView, record::Numeric, Read};

use crate::errors::Error;
use crate::utils::genomics::calculate_dynamic_svlen;
use crate::utils::stats::phred_to_prob;

/* ============ Data Structures =================== */
/// A struct that records samples with their index map.
#[derive(Debug, Clone)]
pub(crate) struct SampleInfo {
    pub samples: Vec<String>,
    pub samples_index_map: HashMap<String, usize>,
}
/* ================================================ */

/* ========= BCF Extraction Functions ============= */
/// Extract sample names from VCF header
pub(crate) fn extract_sample_names(vcf: &bcf::Reader) -> Vec<String> {
    let header = vcf.header();
    header
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect()
}

/// Get chromosome name from a VCF record
pub(crate) fn get_chrom(record: &bcf::Record, header: &HeaderView) -> Result<String> {
    let rid = record
        .rid()
        .ok_or_else(|| Error::VcfRecordMissingChromosome {
            chrom: "unknown".to_string(),
            pos: record.pos(),
        })?;

    let chrom_bytes = header
        .rid2name(rid)
        .map_err(|_| Error::VcfChromResolutionFailed {
            pos: record.pos(),
            rid,
            details: "Failed to resolve chromosome name".to_string(),
        })?;

    let chrom = String::from_utf8_lossy(chrom_bytes).to_string();

    Ok(chrom)
}

/// Get SVLEN from INFO field or calculate dynamically
/// Positive = insertion, Negative = deletion
/// Example: REF=ACAG, ALT=ACAGCAG → SVLEN=+3
pub(crate) fn get_svlen(
    record: &bcf::Record,
    alt_idx: usize,
    ref_seq: &[u8],
    alt_seq: &[u8],
) -> Result<i32> {
    // Try to get SVLEN from INFO field first
    if let Ok(Some(svlens)) = record.info(b"SVLEN").integer() {
        if let Some(&svlen) = svlens.get(alt_idx) {
            if !svlen.is_missing() {
                return Ok(svlen);
            }
        }
    }

    // Fallback: calculate dynamically using anchor detection
    Ok(calculate_dynamic_svlen(ref_seq, alt_seq))
}

// @TODO: Currently caters Phred and Linear Prob.
// Check if LogProb is also possible and extend accordingly.
/// Get combined probability that variant is absent (artifact)
///
/// Extracts and combines two sources of artifact evidence:
/// 1. PROBABILITY_ABSENT tag: Primary absent probability
/// 2. PROBABILITY_ARTIFACT tag: Additional artifact probability
///
/// NOTE: Both can be either PHRED-scaled or direct probabilities.
///
/// # Formula
/// ```text
/// P(absent) = P(artifact_PA) + P(artifact_ART)
/// ```
///
/// # Arguments
/// * `record` - VCF record
/// * `header` - VCF header (for error messages)
/// * `sample_idx` - Sample index (0-based)
/// * `is_phred` - Whether probabilities are PHRED-scaled
///
/// # Returns
/// - `Some(prob)` if both tag exists
/// - `None` if both tags are missing
///
/// # Errors
/// Returns error if:
/// - Probability is NaN
/// - Final combined probability is outside [0, 1]
pub(crate) fn get_prob_absent(
    record: &bcf::Record,
    header: &HeaderView,
    alt_idx: usize,
    is_phred: bool,
) -> Result<Option<f64>> {
    let prob_absent = match record.info(b"PROB_ABSENT").float()? {
        Some(p) if alt_idx < p.len() => p[alt_idx],
        _ => return Ok(None),
    };

    let prob_artifact = match record.info(b"PROB_ARTIFACT").float()? {
        Some(p) if alt_idx < p.len() => p[alt_idx],
        _ => return Ok(None),
    };

    if prob_absent.is_missing() || prob_artifact.is_missing() {
        return Ok(None);
    }

    if prob_absent.is_nan() {
        return Err(Error::InvalidProbabilityValue {
            field: "PROB_ABSENT".to_string(),
            value: prob_absent,
            chrom: get_chrom(record, header)?,
            pos: record.pos(),
        }
        .into());
    }

    if prob_artifact.is_nan() {
        return Err(Error::InvalidProbabilityValue {
            field: "PROB_ARTIFACT".to_string(),
            value: prob_artifact,
            chrom: get_chrom(record, header)?,
            pos: record.pos(),
        }
        .into());
    }

    let probability_absent = if is_phred {
        phred_to_prob(prob_absent as f64)
    } else {
        prob_absent as f64
    };

    let probability_artifact = if is_phred {
        phred_to_prob(prob_artifact as f64)
    } else {
        prob_artifact as f64
    };

    let probability = probability_absent + probability_artifact;

    if !(0.0..=1.0).contains(&probability) {
        return Err(Error::InvalidProbabilityValue {
            field: "DERIVED PROBABILITY ABSENT: PROB_ABSENT + PROBABILITY_ARTIFACT".to_string(),
            value: probability as f32,
            chrom: get_chrom(record, header)?,
            pos: record.pos(),
        }
        .into());
    }

    Ok(Some(probability))
}

/// Extract sample allele frequencies for specific alt allele
pub(crate) fn get_sample_afs(
    record: &bcf::Record,
    header: &HeaderView,
    samples_index_map: &HashMap<String, usize>,
    alt_idx: usize,
) -> Result<HashMap<String, f64>> {
    let mut sample_afs = HashMap::new();

    let afs = match record.format(b"AF").float() {
        Ok(a) => a,
        Err(_) => {
            warn!(
                "AF field missing at {}:{} - variant will have no AF data",
                get_chrom(record, header)?,
                record.pos() + 1
            );
            return Ok(sample_afs);
        }
    };

    for (sample_name, &vcf_header_idx) in samples_index_map {
        let Some(sample_af_values) = afs.get(vcf_header_idx) else {
            continue;
        };

        let Some(&af) = sample_af_values.get(alt_idx) else {
            continue;
        };

        if af.is_missing() {
            continue;
        }

        if af.is_nan() || !(0.0..=1.0).contains(&af) {
            return Err(Error::InvalidAlleleFrequency {
                sample: sample_name.clone(),
                af,
                chrom: get_chrom(record, header)?,
                pos: record.pos(),
            }
            .into());
        }

        sample_afs.insert(sample_name.clone(), af as f64);
    }

    if sample_afs.is_empty() {
        debug!(
            "No valid AF values for any sample at {}:{} - variant will be skipped in analysis",
            get_chrom(record, header)?,
            record.pos() + 1
        );
    }

    Ok(sample_afs)
}
/* ================================================ */

/* ======== BCF Allele Type Check Functions ======= */
// /// Check if all alleles are SNVs (all same length as ref)
// /// Example: REF=A, ALT=T,G -> all length 1
// NOTE: Not required, can be toggled on if required in other utilities.
// pub fn all_alleles_snv(ref_allele: &[u8], alt_alleles: &[&[u8]]) -> bool {
//     let ref_allele_len = ref_allele.len();
//     alt_alleles.iter().all(|alt| alt.len() == ref_allele_len)
// }

/// Check if allele represents no variant (reference)
///
/// Matches alleles that indicate "no alternative allele":
/// - `.` - Missing/no ALT allele (VCF spec)
/// - `<REF>` - Explicit reference allele (rare)
pub(crate) fn is_reference_allele(allele: &[u8]) -> bool {
    allele == b"." || allele == b"<REF>"
}

/// Check if allele is an indel (different length than ref)
/// Example: REF=ACAG, ALT=ACAGCAG -> indel
pub(crate) fn is_indel(ref_allele: &[u8], alt_allele: &[u8]) -> bool {
    ref_allele.len() != alt_allele.len()
}

/// Check if allele is symbolic (starts with <)
/// Examples:
/// - `<DEL>` - Deletion
/// - `<INS>` - Insertion  
/// - `<DUP>` - Duplication
/// - `<INV>` - Inversion
/// - `<CNV>` - Copy number variation
/// - `<NON_REF>` - GVCF placeholder for "any other allele"
pub(crate) fn is_symbolic(allele: &[u8]) -> bool {
    allele.len() >= 3 && allele.starts_with(b"<") && allele.ends_with(b">")
}

/// Check if allele is a breakend (contains [ or ])
/// Example: ALT=A[chr2:100[, ]chr3:200]T
pub(crate) fn is_breakend(allele: &[u8]) -> bool {
    allele.iter().any(|&c| c == b'[' || c == b']')
}

/// Check if allele is a spanning deletion (*)
/// Example: ALT=*
pub(crate) fn is_spanning_deletion(allele: &[u8]) -> bool {
    allele == b"*"
}
/* ================================================ */

/* ========= BCF Validation Functions ============= */
/// Validate VCF file format by checking samples and first record
/// and returning final samples list to process after exclusion of
/// those that user wanted to get excluded from processing.
pub(crate) fn validate_vcf_file(
    vcf_path: &PathBuf,
    samples_exclusion: &[String],
) -> Result<SampleInfo> {
    info!("Validating VCF file format: {}", vcf_path.display());

    let mut vcf = bcf::Reader::from_path(vcf_path).context("Failed to open VCF/BCF file")?;

    let sample_names = extract_sample_names(&vcf);

    if sample_names.is_empty() {
        return Err(Error::VcfNoSamples.into());
    }

    let mut invalid_exclusions: Vec<String> = Vec::new();

    for excluded_sample in samples_exclusion {
        if !sample_names.contains(excluded_sample) {
            invalid_exclusions.push(excluded_sample.clone());
        }
    }

    if !invalid_exclusions.is_empty() {
        return Err(Error::InvalidSampleExclusion {
            samples: invalid_exclusions.join(", "),
        }
        .into());
    }

    let mut remaining_samples: Vec<String> = vec![];
    let mut samples_index_map: HashMap<String, usize> = HashMap::new();

    for (i, s) in sample_names.iter().enumerate() {
        if !samples_exclusion.contains(s) {
            remaining_samples.push(s.to_string());
            samples_index_map.insert(s.to_string(), i);
        }
    }

    if remaining_samples.is_empty() {
        return Err(Error::NoSamplesAfterExclusion.into());
    }

    info!("Samples to process: {}", remaining_samples.len());
    info!("Sample names: {:?}", remaining_samples);

    let header = vcf.header().clone();

    match vcf.records().next() {
        None => {
            return Err(Error::VcfEmpty.into());
        }
        Some(Err(e)) => {
            return Err(Error::VcfRecordRead {
                details: e.to_string(),
            }
            .into());
        }
        Some(Ok(record)) => {
            let chrom = get_chrom(&record, &header)?;
            let pos = record.pos();
            info!("  First variant: {}:{}", chrom, pos + 1);
        }
    }

    info!("VCF file format validated successfully");

    Ok(SampleInfo {
        samples: remaining_samples,
        samples_index_map,
    })
}

/* ================================================ */

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;
    use crate::utils::stats::TEST_EPSILON;


    #[test]
    fn test_is_reference_allele() {
        assert!(is_reference_allele(b"."));
        assert!(is_reference_allele(b"<REF>"));
    }

    #[test]
    fn test_is_indel() {
        assert!(is_indel(b"AC", b"ACG")); // insertion
        assert!(is_indel(b"ACG", b"AC")); // deletion
        assert!(!is_indel(b"AC", b"AG")); // SNV
    }

    #[test]
    fn test_is_symbolic() {
        // Structural variants
        assert!(is_symbolic(b"<DEL>"));
        assert!(is_symbolic(b"<INS>"));
        assert!(is_symbolic(b"<DUP>"));
        assert!(is_symbolic(b"<INV>"));
        assert!(is_symbolic(b"<CNV>"));
        assert!(is_symbolic(b"<DUP:TANDEM>"));

        // GVCF placeholder
        assert!(is_symbolic(b"<NON_REF>"));
    }

    #[test]
    fn test_is_breakend() {
        assert!(is_breakend(b"A[chr2:100["));
        assert!(is_breakend(b"]chr3:200]T"));
        assert!(!is_breakend(b"ACG"));
    }

    #[test]
    fn test_is_spanning_deletion() {
        assert!(is_spanning_deletion(b"*"));
        assert!(!is_spanning_deletion(b"AC"));
    }

    fn encode_genotype_allele(allele_index: i32, phased: bool) -> i32 {
        let phased_flag = if phased { 1 } else { 0 };
        (allele_index + 1) * 2 | phased_flag
    }

    fn create_test_vcf() -> (NamedTempFile, Vec<String>) {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        let mut header = rust_htslib::bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(br##"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">"##);
        header.push_record(
            br##"##INFO=<ID=PROB_ABSENT,Number=1,Type=Float,Description="Probability absent">"##,
        );
        header.push_record(
            br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Probability artifact">"##,
        );
        header.push_sample(b"sample1");
        header.push_sample(b"sample2");
        header.push_record(br##"##FORMAT=<ID=GT,Number=.,Type=String,Description="Genotype">"##);
        header.push_record(
            br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">"##,
        );

        let mut wtr = rust_htslib::bcf::Writer::from_path(
            path,
            &header,
            false,
            rust_htslib::bcf::Format::Vcf,
        )
        .unwrap();

        let mut rec = wtr.empty_record();
        rec.set_rid(Some(0)); // RID 0 corresponds to contig chr1
        rec.set_pos(99); // 0-based position 99 to 100 in VCF
        rec.set_alleles(&[b"A", b"T"]).unwrap();

        rec.push_info_integer(b"SVLEN", &[1]).unwrap();
        rec.push_info_float(b"PROB_ABSENT", &[0.01]).unwrap();
        rec.push_info_float(b"PROB_ARTIFACT", &[0.005]).unwrap();

        let genotypes_data = [
            encode_genotype_allele(0, false),
            encode_genotype_allele(1, false),
            encode_genotype_allele(1, false),
            encode_genotype_allele(1, false),
        ];

        rec.push_format_integer(b"GT", &genotypes_data).unwrap();

        let af_values = [0.45_f32, 0.98_f32];
        rec.push_format_float(b"AF", &af_values).unwrap();

        wtr.write(&rec).unwrap();

        (tmp, vec!["sample1".to_string(), "sample2".to_string()])
    }

    #[test]
    fn test_validate_vcf_file() {
        let (tmp_vcf, _) = create_test_vcf();

        // Test with no exclusions
        let samples_info = validate_vcf_file(&tmp_vcf.path().to_path_buf(), &[]).unwrap();
        assert_eq!(samples_info.samples.len(), 2);
        assert_eq!(
            samples_info.samples,
            vec!["sample1".to_string(), "sample2".to_string()]
        );
        assert_eq!(samples_info.samples_index_map["sample1"], 0);
        assert_eq!(samples_info.samples_index_map["sample2"], 1);

        // Test with one exclusion
        let samples_info =
            validate_vcf_file(&tmp_vcf.path().to_path_buf(), &["sample1".to_string()]).unwrap();
        assert_eq!(samples_info.samples, vec!["sample2".to_string()]);
        assert_eq!(samples_info.samples_index_map["sample2"], 1);

        // // Test exclusion of non-existent sample
        let err = validate_vcf_file(
            &tmp_vcf.path().to_path_buf(),
            &["invalid_sample".to_string()],
        )
        .unwrap_err();
        assert!(format!("{err:?}").contains("invalid_sample"));
    }

    #[test]
    fn test_get_chrom_and_svlen() {
        let (tmp_vcf, _) = create_test_vcf();
        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let chrom = get_chrom(&record, &header).unwrap();
        assert_eq!(chrom, "chr1"); // because rid=0 corresponds to first header sample

        let svlen = get_svlen(&record, 0, b"A", b"T").unwrap();
        assert_eq!(svlen, 1);
    }

    #[test]
    fn test_get_prob_absent() {
        let (tmp_vcf, _) = create_test_vcf();
        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let prob = get_prob_absent(&record, &header, 0, false)
            .unwrap()
            .unwrap();
        assert!((prob - 0.015).abs() < 1e-6);
    }

    //add test_get_prob_phred
    // test_get_prob_absent_missing

    #[test]
    fn test_get_sample_afs_all_samples() {
        let (tmp_vcf, _) = create_test_vcf();

        // Process all samples
        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);
        samples_index_map.insert("sample2".to_string(), 1);

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let sample_afs = get_sample_afs(&record, &header, &samples_index_map, 0).unwrap();

        // Should have both samples
        assert_eq!(sample_afs.len(), 2);
        assert!((sample_afs["sample1"] - 0.45).abs() < 1e-6);
        assert!((sample_afs["sample2"] - 0.98).abs() < 1e-6);
    }

    #[test]
    fn test_get_sample_afs_all_missing() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        let mut header = rust_htslib::bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_sample(b"sample1");
        header.push_sample(b"sample2");
        header.push_record(br##"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"##);
        header.push_record(
            br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">"##,
        );

        let mut wtr = rust_htslib::bcf::Writer::from_path(
            path,
            &header,
            false,
            rust_htslib::bcf::Format::Vcf,
        )
        .unwrap();

        let mut rec = wtr.empty_record();
        rec.set_rid(Some(0));
        rec.set_pos(99);
        rec.set_alleles(&[b"ACAG", b"ACAGCAG"]).unwrap();

        let missing = f32::missing();
        let af_values = [missing, missing]; // Both samples missing
        rec.push_format_float(b"AF", &af_values).unwrap();

        wtr.write(&rec).unwrap();
        drop(wtr);

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);
        samples_index_map.insert("sample2".to_string(), 1);

        let mut reader = bcf::Reader::from_path(path).unwrap();
        let header_view = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let sample_afs = get_sample_afs(&record, &header_view, &samples_index_map, 0).unwrap();

        assert!(sample_afs.is_empty());
        assert_eq!(sample_afs.len(), 0);
    }

    #[test]
    fn test_get_sample_afs_with_exclusion() {
        let (tmp_vcf, _) = create_test_vcf();

        // Exclude sample1 - only process sample2
        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample2".to_string(), 1);

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let sample_afs = get_sample_afs(&record, &header, &samples_index_map, 0).unwrap();

        // Should only have sample2, NOT sample1
        assert_eq!(sample_afs.len(), 1);
        assert!(!sample_afs.contains_key("sample1"));
        assert!(sample_afs.contains_key("sample2"));
        assert!((sample_afs["sample2"] - 0.98).abs() < 1e-6);
    }
}

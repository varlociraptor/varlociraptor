//! bcf_utils.rs
//!
//! Utilities for VCF/BCF file handling and parsing.
//!
//! This module provides:
//! 1. Sample information extraction
//! 2. Record field extraction (chromosome, SVLEN, probabilities, allele frequencies)
//! 3. Allele type classification (indel, symbolic, breakend, reference, spanning deletion)
//! 4. VCF file validation

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, Result};
use log::{info, warn};
use rust_htslib::bcf::{self, header::HeaderView, record::Numeric, Read};

use crate::errors::Error;
use crate::utils::genomics::calculate_dynamic_svlen;
use crate::utils::stats::phred_to_prob;

/* ============ Data Structures =================== */

/// Sample information extracted from VCF header.
///
/// Contains sample names and their corresponding indices
/// in the VCF header for efficient lookup during processing.
#[derive(Debug, Clone)]
pub(crate) struct SampleInfo {
    /// Sample names to process
    pub samples: Vec<String>,
    /// Map of sample name to VCF header index
    pub samples_index_map: HashMap<String, usize>,
}

/* ================================================ */

/* ========= BCF Extraction Functions ============= */

/// Extract sample names from VCF header
///
/// # Arguments
/// * `vcf` - VCF reader
///
/// # Returns
/// Vector of sample names as Strings
///
/// # Example
/// assert_eq!(extract_sample_names(&vcf), vec!["sample1", "sample2"]);
pub(crate) fn extract_sample_names(vcf: &bcf::Reader) -> Vec<String> {
    let header = vcf.header();
    header
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect()
}

/// Get chromosome name from a VCF record
///
/// # Arguments
/// * `record` - VCF record
/// * `header` - VCF header (for resolving RID to name)
///
/// # Returns
/// Chromosome name as String
///
/// # Errors
/// Returns error if:
/// - RID is missing in record
/// - Chromosome name resolution fails
///
/// # Example
/// assert_eq!(get_chrom(&record, &header).unwrap(), "chr1");
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

/// Get SVLEN from INFO field or calculate dynamically.
///
/// Returns length difference between ALT and REF alleles.
/// (above referes to SVLEN, not simple length difference)
/// Positive = insertion, Negative = deletion.
///
/// # Arguments
/// * `record` - VCF record
/// * `alt_idx` - Index into ALT alleles (0 = first ALT)
/// * `ref_seq` - Reference allele sequence
/// * `alt_seq` - Alternate allele sequence
///
/// # Example
/// REF=ACAG, ALT=ACAGCAG : SVLEN=+3
/// assert_eq!(get_svlen(&record, 0, b"ACAG", b"ACAGCAG").unwrap(), 3);
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
/// Notes: Currently caters Phred and Linear Prob.
/// Check if LogProb is also possible and extend accordingly if required.
///
/// # Arguments
/// * `record` - VCF record
/// * `header` - VCF header (for error messages)
/// * `alt_idx` - ALT allele index (0 = first ALT)
/// * `is_phred` - Whether probabilities are PHRED-scaled
///
/// # Returns
/// - `Some(prob)` if both PROB_ABSENT and PROB_ARTIFACT exist for this ALT
/// - `None` if either tag is missing or missing at this ALT index
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

/// Extract per-sample allele frequencies for a specific ALT allele.
///
/// Reads FORMAT:AF field and returns AF values for each sample
/// in the samples_index_map.
///
/// # Arguments
/// * `record` - VCF record
/// * `header` - VCF header
/// * `samples_index_map` - Map of sample names to VCF indices
/// * `alt_idx` - Index into ALT alleles (0 = first ALT)
///
/// # Returns
/// HashMap of sample name to AF value. Empty if AF field missing.
///
/// # Errors
/// Returns error if AF value is NaN or outside [0.0, 1.0].
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
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is reference, `false` otherwise
///
/// # Example
/// assert!(is_reference_allele(b"."));
pub(crate) fn is_reference_allele(allele: &[u8]) -> bool {
    allele == b"." || allele == b"<REF>"
}

/// Check if allele is an indel (different length than ref)
///
/// # Arguments
/// * `ref_allele` - Reference allele sequence
/// * `alt_allele` - Alternate allele sequence
///
/// # Returns
/// `true` if allele is an indel, `false` otherwise
///
/// # Example
/// assert!(is_indel(b"ACAG", b"ACAGCAG"));
pub(crate) fn is_indel(ref_allele: &[u8], alt_allele: &[u8]) -> bool {
    ref_allele.len() != alt_allele.len()
}

/// Check if allele is symbolic (starts with <)
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is symbolic, `false` otherwise
///
/// Examples:
/// - `<DEL>` - Deletion
/// - `<INS>` - Insertion  
/// - `<DUP>` - Duplication
/// - `<INV>` - Inversion
/// - `<CNV>` - Copy number variation
/// - `<NON_REF>` - GVCF placeholder for "any other allele"
/// assert!(is_symbolic(b"<DEL>"));
pub(crate) fn is_symbolic(allele: &[u8]) -> bool {
    allele.len() >= 3 && allele.starts_with(b"<") && allele.ends_with(b">")
}

/// Check if allele is a breakend (contains [ or ])
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is a breakend, `false` otherwise
///
/// # Examples
/// assert!(is_breakend(b"A[chr2:100["));
pub(crate) fn is_breakend(allele: &[u8]) -> bool {
    allele.iter().any(|&c| c == b'[' || c == b']')
}

/// Check if allele is a spanning deletion (*)
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is a spanning deletion, `false` otherwise
///
/// # Example
/// assert!(is_spanning_deletion(b"*"));
pub(crate) fn is_spanning_deletion(allele: &[u8]) -> bool {
    allele == b"*"
}

/* ================================================ */

/* ========= BCF Validation Functions ============= */

/// Validate VCF file and extract sample information.
///
/// Performs validation checks:
/// 1. File can be opened
/// 2. Contains at least one sample
/// 3. Excluded samples exist in file
/// 4. At least one sample remains after exclusion
/// 5. Contains at least one variant record
///
/// # Arguments
/// * `vcf_path` - Path to VCF/BCF file
/// * `samples_exclusion` - Sample names to exclude from processing
///
/// # Returns
/// `SampleInfo` with remaining samples and their indices
///
/// # Example
/// assert!(validate_vcf_file(&vcf_path, &vec!["sample1".to_string()]).is_ok());
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
    use rust_htslib::bcf::record::Numeric;
    use tempfile::NamedTempFile;

    use crate::utils::stats::TEST_EPSILON;

    /// Encodes a single allele for the GT field in BCF format.
    ///
    /// BCF stores each allele as: (allele_index + 1) << 1 | phased_bit
    /// - allele_index: -1 for missing, 0 for REF, 1 for ALT1, etc.
    /// - phased_bit: 1 for phased ('|'), 0 for unphased ('/')
    ///
    /// Example:
    ///   encode_genotype_allele(0, false) -> 2   // allele 0, unphased
    ///
    /// A full genotype like "0/0" would be encoded as: [2, 2]
    ///
    /// Reference: https://samtools.github.io/hts-specs/BCFv2_qref.pdf
    fn encode_genotype_allele(allele_index: i32, phased: bool) -> i32 {
        let phased_flag = if phased { 1 } else { 0 };
        (allele_index + 1) * 2 | phased_flag
    }

    /// Configuration for test VCF creation
    struct TestVcfConfig<'a> {
        ref_allele: &'a [u8],
        alt_alleles: Vec<&'a [u8]>,
        af_values: Option<Vec<f32>>,
        prob_absent: Option<Vec<f32>>,
        prob_artifact: Option<Vec<f32>>,
        num_samples: usize,
        use_phred: bool,
    }

    impl<'a> Default for TestVcfConfig<'a> {
        fn default() -> Self {
            Self {
                ref_allele: b"A",
                alt_alleles: vec![b"AT"],
                af_values: None,
                prob_absent: None,
                prob_artifact: None,
                num_samples: 2,
                use_phred: false,
            }
        }
    }

    fn create_test_vcf(config: TestVcfConfig) -> (NamedTempFile, Vec<String>) {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        let mut header = rust_htslib::bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##);

        // Conditionally set PROB_ABSENT and PROB_ARTIFACT headers based on use_phred
        if config.use_phred {
            header.push_record(br##"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,Description="Probability absent (PHRED)">"##);
            header.push_record(br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Probability artifact (PHRED)">"##);
        } else {
            header.push_record(br##"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,Description="Probability absent (linear)">"##);
            header.push_record(br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Probability artifact (linear)">"##);
        }

        header.push_record(br##"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"##);
        header.push_record(
            br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">"##,
        );

        let mut sample_names = Vec::new();
        for i in 0..config.num_samples {
            let name = format!("sample{}", i + 1);
            header.push_sample(name.as_bytes());
            sample_names.push(name);
        }

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

        // Build alleles: [REF, ALT1, ALT2, ...]
        let mut alleles_vec = vec![config.ref_allele];
        alleles_vec.extend(config.alt_alleles.iter());
        rec.set_alleles(&alleles_vec).unwrap();

        // SVLEN for each ALT
        let svlen_values: Vec<i32> = config
            .alt_alleles
            .iter()
            .map(|alt| (alt.len() as i32) - (config.ref_allele.len() as i32))
            .collect();
        rec.push_info_integer(b"SVLEN", &svlen_values).unwrap();

        // PROB_ABSENT
        let prob_absent = match config.prob_absent {
            Some(pa) => pa,
            None => {
                if config.use_phred {
                    vec![20.0; config.alt_alleles.len()]
                } else {
                    vec![0.005; config.alt_alleles.len()]
                }
            }
        };
        rec.push_info_float(b"PROB_ABSENT", &prob_absent).unwrap();

        // PROB_ARTIFACT
        let prob_artifact = match config.prob_artifact {
            Some(pa) => pa,
            None => {
                if config.use_phred {
                    vec![10.0; config.alt_alleles.len()]
                } else {
                    vec![0.005; config.alt_alleles.len()]
                }
            }
        };
        rec.push_info_float(b"PROB_ARTIFACT", &prob_artifact)
            .unwrap();

        // GT
        let mut genotypes = Vec::new();
        for _ in 0..config.num_samples {
            genotypes.push(encode_genotype_allele(0, false));
            genotypes.push(encode_genotype_allele(1, false));
        }
        rec.push_format_integer(b"GT", &genotypes).unwrap();

        // AF values
        if let Some(afs) = config.af_values {
            rec.push_format_float(b"AF", &afs).unwrap();
        }

        wtr.write(&rec).unwrap();

        (tmp, sample_names)
    }

    /* ==== BCF Extraction Function(s) tests ========= */

    #[test]
    fn test_get_chrom_and_svlen() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig::default());

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let chrom = get_chrom(&record, &header).unwrap();
        assert_eq!(chrom, "chr1");

        let svlen = get_svlen(&record, 0, b"A", b"AT").unwrap();
        assert_eq!(svlen, 1); // A -> AT is an insertion, svlen = 1
    }

    #[test]
    fn test_get_prob_absent() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            prob_absent: Some(vec![0.01]),
            prob_artifact: Some(vec![0.005]),
            num_samples: 1,
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let prob = get_prob_absent(&record, &header, 0, false)
            .unwrap()
            .unwrap();
        assert!((prob - 0.015).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_get_prob_absent_phred() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            use_phred: true,
            num_samples: 1,
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let prob = get_prob_absent(&record, &header, 0, true).unwrap().unwrap();
        // PHRED 20 ≈ 0.01, PHRED 10 ≈ 0.1, sum ≈ 0.11
        assert!((prob - 0.11).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_get_prob_absent_invalid_sum() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            prob_absent: Some(vec![0.7]),
            prob_artifact: Some(vec![0.6]), // Sum > 1.0
            num_samples: 1,
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let result = get_prob_absent(&record, &header, 0, false);
        assert!(result.is_err());
    }

    #[test]
    fn test_get_sample_afs() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![0.45, 0.98]),
            ..Default::default()
        });

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);
        samples_index_map.insert("sample2".to_string(), 1);

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let sample_afs = get_sample_afs(&record, &header, &samples_index_map, 0).unwrap();

        assert_eq!(sample_afs.len(), 2);
        assert!((sample_afs["sample1"] - 0.45).abs() < TEST_EPSILON);
        assert!((sample_afs["sample2"] - 0.98).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_get_sample_afs_invalid() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![1.5]), // Invalid: > 1.0
            num_samples: 1,
            ..Default::default()
        });

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let result = get_sample_afs(&record, &header, &samples_index_map, 0);
        assert!(result.is_err());
    }

    #[test]
    fn test_get_sample_afs_missing() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![f32::missing()]),
            num_samples: 1,
            ..Default::default()
        });

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let sample_afs = get_sample_afs(&record, &header, &samples_index_map, 0).unwrap();
        assert!(sample_afs.is_empty());
    }

    /* ====== BCF ALLELE Type tests ================== */

    #[test]
    fn test_allele_type_checks() {
        // Reference allele
        assert!(is_reference_allele(b"."));
        assert!(is_reference_allele(b"<REF>"));

        // Indel
        assert!(is_indel(b"AC", b"ACG")); // insertion
        assert!(is_indel(b"ACG", b"AC")); // deletion
        assert!(!is_indel(b"AC", b"AG"));

        // Symbolic
        assert!(is_symbolic(b"<DEL>"));
        assert!(is_symbolic(b"<DUP:TANDEM>"));
        assert!(!is_symbolic(b"ACG"));

        // Breakend
        assert!(is_breakend(b"A[chr2:100["));
        assert!(!is_breakend(b"ACG"));

        // Spanning deletion
        assert!(is_spanning_deletion(b"*"));
        assert!(!is_spanning_deletion(b"AC"));
    }

    /* ======== validate_vcf_file tests ============== */

    #[test]
    fn test_validate_vcf_file() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![0.45, 0.98]),
            ..Default::default()
        });

        // Test with no exclusions
        let samples_info = validate_vcf_file(&tmp_vcf.path().to_path_buf(), &[]).unwrap();
        assert_eq!(samples_info.samples.len(), 2);
        assert_eq!(samples_info.samples_index_map["sample1"], 0);
        assert_eq!(samples_info.samples_index_map["sample2"], 1);

        // Test with one exclusion
        let samples_info =
            validate_vcf_file(&tmp_vcf.path().to_path_buf(), &["sample1".to_string()]).unwrap();
        assert_eq!(samples_info.samples, vec!["sample2".to_string()]);

        // Test exclusion of non-existent sample
        let err = validate_vcf_file(
            &tmp_vcf.path().to_path_buf(),
            &["invalid_sample".to_string()],
        )
        .unwrap_err();
        assert!(format!("{err:?}").contains("invalid_sample"));
    }
}

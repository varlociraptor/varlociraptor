//! mod.rs
//!
//! Microsatellite instability (MSI) preprocessing module.
//! Enriches candidate VCF with:
//! - Dummy indel records for MS regions without a perfect indel
//! - MS region annotations (INFO/REGION_ID)
//!
//! **Note: This feature is experimental.**
//!
//! Owns `PreprocessMSIConfig` (CLI configuration) and
//! `preprocess_ms_candidates` (the entry point that orchestrates header,
//! intersection, variant analysis, writer). Declares these child submodules:
//! 1. `header` - builds the output header: copies INFO field
//!    declarations, adds REGION_ID/MSI_DUMMY
//! 2. `intersection` - streams BED regions against VCF variants,
//!    matching variants to the regions they overlap and adding dummy
//!    indels for regions without a perfect indel
//! 3. `writer` - writes unannotated/annotated real variants and the
//!    dummy deletion and insertion records to the output file
//! 4. `variant_analysis` - classifies whether an indel is a perfect
//!    tandem repeat of the BED region's motif

mod header;
mod intersection;
mod variant_analysis;
mod writer;

use std::path::PathBuf;

use anyhow::{Context, Result};
use rust_htslib::bcf::{self, Format, Read};

use crate::cli::{validate_bed_file, validate_vcf_file};
use crate::utils::aux_info::AuxInfoCollector;
use crate::utils::bcf_utils::{
    validate_info_fields_exist, validate_vcf_file as validate_ms_vcf_file,
};
use crate::utils::ms_bed::validate_bed_file as validate_ms_bed_file;

use header::prepare_header;
use intersection::process_and_annotate;

/* ============ Data Structures =================== */

/// Configuration for MS candidate preprocessing pipeline.
///
/// Contains all parameters needed for the MSI preprocessing workflow,
/// populated from CLI arguments.
pub struct PreprocessMSIConfig {
    /// Path to BED file (sorted) with microsatellite regions.
    pub microsatellite_bed: PathBuf,
    /// Path to candidate VCF/BCF file (sorted) with variant calls (gnomAD-annotated with INFO/HETEROZYGOSITY).
    pub candidate_vcf: PathBuf,
    /// Output file path. Format/compression is inferred from the extension
    /// (see `determine_output_format`'s doc for the full list). If omitted,
    /// writes uncompressed VCF to STDOUT.
    /// Note: output BCF/VCF.GZ is always compressed; there is currently no
    /// option to request uncompressed BCF.
    pub output: Option<PathBuf>,
    /// List of INFO fields to propagate from the candidate VCF to the output.
    pub propagate_info_fields: Vec<String>,
}

impl PreprocessMSIConfig {
    /// Validate the configuration before preprocessing starts.
    ///
    /// Checks run in this order and stop at the first failure: BED extension,
    /// VCF extension, VCF opens, BED first record, VCF first record, requested
    /// INFO fields declared in the VCF header.
    ///
    /// # Returns
    /// * `Ok(())` if every check passes
    /// * `Err(Error::BedFileInvalid)` if the BED path has no `.bed` extension
    /// * `Err(Error::VcfFileInvalid)` if the VCF path has no `.vcf`, `.vcf.gz`, `.vcf.bgz`, `.bcf` or `.bcf.gz` suffix
    /// * `Err` with context `Failed to open candidate VCF` if the VCF cannot be opened
    /// * `Err` from `ms_bed::validate_bed_file` if the BED content check fails
    /// * `Err` from `bcf_utils::validate_vcf_file` if the VCF content check fails
    /// * `Err(Error::VcfInfoFieldsMissing)` if a field in `propagate_info_fields` is not declared in the VCF header
    pub(crate) fn validate(&self) -> Result<()> {
        // Use existing validators from cli.rs to validate file formats.
        validate_bed_file(&self.microsatellite_bed)?;
        validate_vcf_file(&self.candidate_vcf)?;

        let mut vcf =
            bcf::Reader::from_path(&self.candidate_vcf).context("Failed to open candidate VCF")?;

        // Additional file specific validators:
        validate_ms_bed_file(&self.microsatellite_bed)?;
        info!("Validating VCF file: {}", self.candidate_vcf.display());
        validate_ms_vcf_file(&mut vcf)?;
        validate_info_fields_exist(vcf.header(), &self.propagate_info_fields)?;

        Ok(())
    }

    /// Determine output format and compression from file extension.
    ///
    /// # Format Detection
    /// - `.bcf` to BCF binary, compressed
    /// - `.bcf.gz` to BCF binary, compressed (identical to `.bcf`)
    /// - `.vcf.gz` or `.vcf.bgz` to VCF text, BGZF compressed
    /// - `.vcf` to VCF text, uncompressed
    /// - Anything else to VCF text, uncompressed (safest default)
    /// - None (stdout) to VCF text, uncompressed
    ///
    /// # Returns
    /// Tuple of (`Format`, `uncompressed`) for writer creation, `true` means no compression is applied
    fn determine_output_format(&self) -> (Format, bool) {
        match self.output.as_deref() {
            None => (Format::Vcf, true),
            Some(path) => {
                let path_str = path.to_string_lossy().to_lowercase();

                if path_str.ends_with(".bcf") || path_str.ends_with(".bcf.gz") {
                    // .bcf.gz produces identical output to .bcf.
                    (Format::Bcf, false)
                } else if path_str.ends_with(".vcf.gz") || path_str.ends_with(".vcf.bgz") {
                    (Format::Vcf, false)
                } else {
                    // .vcf and any other extensions default to uncompressed VCF
                    (Format::Vcf, true)
                }
            }
        }
    }

    /// Log the input paths, the output destination and the INFO fields to propagate.
    ///
    /// # Returns
    /// Nothing, as all output goes through `info!`
    fn log_config(&self) {
        info!("Input files:");
        info!("BED file: {}", self.microsatellite_bed.display());
        info!("VCF/BCF file: {}", self.candidate_vcf.display());

        if let Some(ref output) = self.output {
            info!("Output: {}", output.display());
        } else {
            info!("Output: STDOUT (VCF format)");
        }

        if self.propagate_info_fields.is_empty() {
            info!(
                "Propagate INFO fields: none (standard fields only, if present: {})",
                msi_omit_aux_info!()
            );
        } else {
            info!(
                "Propagate INFO fields: {}",
                self.propagate_info_fields.join(", ")
            );
        }
    }
}

/* ============ Functions ========================= */

/// Run the MSI candidate preprocessing workflow end to end.
///
/// Opens the candidate VCF, builds the output header, creates the writer
/// (STDOUT if no output path is set), streams the BED regions against the
/// variants via `process_and_annotate`, then logs the final stats.
/// Callers run `PreprocessMSIConfig::validate` first.
///
/// # Arguments
/// * `config` - Preprocessing configuration
///
/// # Returns
/// * `Ok(())` once the output has been written
/// * `Err` with context `Failed to open input VCF` if the candidate VCF cannot be opened
/// * `Err` if `AuxInfoCollector::new` fails
/// * `Err` if `prepare_header` fails
/// * `Err` if the output writer cannot be created
/// * `Err` if `process_and_annotate` fails
pub fn preprocess_ms_candidates(config: PreprocessMSIConfig) -> Result<()> {
    info!("----------------------------------------------");
    info!("Step 1: Configuration");
    info!("----------------------------------------------");
    config.log_config();

    info!("----------------------------------------------");
    info!("Step 2: Streaming Intersection & Output Generation");
    info!("----------------------------------------------");
    info!("Starting step 2.");

    info!("Opening input VCF: {}", config.candidate_vcf.display());
    let mut input_vcf =
        bcf::Reader::from_path(&config.candidate_vcf).context("Failed to open input VCF")?;

    info!("Preparing auxiliary INFO field collector...");
    let aux_fields: Vec<Vec<u8>> = config
        .propagate_info_fields
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
    let aux_info_collector = AuxInfoCollector::new(&aux_fields, &input_vcf)?;

    info!("Preparing output header...");
    let output_header = prepare_header(
        input_vcf.header(),
        &config.microsatellite_bed,
        &aux_info_collector,
    )?;

    info!("Determining output format based on output path...");
    let (output_format, output_uncompressed) = config.determine_output_format();

    info!("Creating output writer...");
    let mut writer = match config.output {
        Some(ref path) => {
            bcf::Writer::from_path(path, &output_header, output_uncompressed, output_format)?
        }
        None => bcf::Writer::from_stdout(&output_header, output_uncompressed, output_format)?,
    };

    info!("\nStarting streaming intersection...");
    let stats = process_and_annotate(
        &mut input_vcf,
        &config.microsatellite_bed,
        &mut writer,
        &aux_info_collector,
    )?;
    info!("Intersection finished and output generated.");

    info!("----------------------------------------------");
    info!("Step 3(Final): Logging Final Stats");
    info!("----------------------------------------------");
    stats.log_stats();
    info!("==============================================");
    info!("Preprocessing MSI complete");
    info!("==============================================");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::errors::Error;

    // Helper: Build a PreprocessMSIConfig with the given output path and dummy inputs
    fn create_test_config(output: Option<PathBuf>) -> PreprocessMSIConfig {
        PreprocessMSIConfig {
            microsatellite_bed: PathBuf::from("test.bed"),
            candidate_vcf: PathBuf::from("test.vcf"),
            output,
            propagate_info_fields: vec![],
        }
    }

    /* ====== determine_output_format tests ========== */

    #[test]
    fn test_determine_output_format_stdout() {
        let config = create_test_config(None);
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(uncompressed, "Stdout should be uncompressed");
    }

    #[test]
    fn test_determine_output_format_bcf() {
        let config = create_test_config(Some(PathBuf::from("output.bcf")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Bcf);
        assert!(!uncompressed, "BCF should be compressed");
    }

    #[test]
    fn test_determine_output_format_bcf_gz() {
        let config = create_test_config(Some(PathBuf::from("output.bcf.gz")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Bcf);
        assert!(!uncompressed, ".bcf.gz should be compressed");
    }

    #[test]
    fn test_determine_output_format_vcf_uncompressed() {
        let config = create_test_config(Some(PathBuf::from("output.vcf")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(uncompressed, "Plain .vcf should be uncompressed");
    }

    #[test]
    fn test_determine_output_format_vcf_gz() {
        let config = create_test_config(Some(PathBuf::from("output.vcf.gz")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(!uncompressed, ".vcf.gz should be BGZF compressed");
    }

    #[test]
    fn test_determine_output_format_vcf_bgz() {
        let config = create_test_config(Some(PathBuf::from("output.vcf.bgz")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(!uncompressed, ".vcf.bgz should be BGZF compressed");
    }

    #[test]
    fn test_determine_output_format_case_insensitive() {
        let config = create_test_config(Some(PathBuf::from("OUTPUT.BCF")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Bcf);
        assert!(!uncompressed);

        let config = create_test_config(Some(PathBuf::from("output.VCF.GZ")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(!uncompressed);
    }

    #[test]
    fn test_determine_output_format_unknown_extension() {
        let config = create_test_config(Some(PathBuf::from("weird.txt")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(
            uncompressed,
            "Unknown extension should default to uncompressed"
        );

        let config = create_test_config(Some(PathBuf::from("output")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(uncompressed, "No extension should default to uncompressed");
    }

    /* ====== validate tests ========================= */

    #[test]
    fn test_validate_bed_wrong_extension() {
        let bed = PathBuf::from("regions.txt");
        let config = PreprocessMSIConfig {
            microsatellite_bed: bed.clone(),
            ..create_test_config(None)
        };

        let err = config.validate().unwrap_err();
        assert_eq!(
            err.downcast_ref::<Error>(),
            Some(&Error::BedFileInvalid { path: bed })
        );
    }

    #[test]
    fn test_validate_vcf_wrong_extension() {
        let vcf = PathBuf::from("calls.txt");
        let config = PreprocessMSIConfig {
            candidate_vcf: vcf.clone(),
            ..create_test_config(None)
        };

        let err = config.validate().unwrap_err();
        assert_eq!(
            err.downcast_ref::<Error>(),
            Some(&Error::VcfFileInvalid { path: vcf })
        );
    }

    #[test]
    fn test_validate_vcf_missing_file() {
        let err = create_test_config(None).validate().unwrap_err();
        assert_eq!(err.to_string(), "Failed to open candidate VCF");
    }
}

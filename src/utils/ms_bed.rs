//! ms_bed.rs
//!
//! BED file parsing utilities for microsatellite instability analysis.
//!
//! This module provides utilities for:
//! 1. Parsing BED records with microsatellite motif information
//! 2. Validating BED file format and content
//!
//! Expected BED format: chrom, start, end, name
//! where name follows the pattern `<count>x<motif>` (e.g., "15xCAG")

use std::path::PathBuf;

use anyhow::{Context, Result};
use bio::io::bed;
use log::info;

use crate::errors::Error;

/// Microsatellite region from BED file.
///
/// Represents a genomic region containing a microsatellite repeat,
/// parsed from BED format with motif information in the name field.
#[derive(Debug, Clone)]
pub(crate) struct BedRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub motif: String,
}

impl BedRegion {
    /// Returns formatted region identifier.
    ///
    /// # Example
    /// "chr1:100-200"
    pub(crate) fn region_id(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start, self.end)
    }

    /// Returns the motif length in bases.
    pub(crate) fn motif_length(&self) -> usize {
        self.motif.len()
    }

    /// Returns true if motif length is valid for MSI analysis (1-6 bases).
    pub(crate) fn is_valid_motif(&self) -> bool {
        (1..=6).contains(&self.motif_length())
    }

    /*  NOTE: flexible is_valid_motif can be implemented if later need
        for MS regions to be greater than length 6.
    */
}

/// Parse motif from BED name field.
///
/// Extracts the motif sequence from format `<count>x<motif>`.
///
/// # Arguments
/// * `name` - BED record name field (e.g., "5xCAG")
///
/// # Returns
/// * `Ok(motif)` - Extracted motif string (e.g., "CAG")
/// * `Err` - Invalid format
///
/// # Errors
/// Returns error if:
/// - Missing 'x' separator
/// - Non-numeric repeat count
/// - Zero repeat count
/// - Empty motif
///
/// # Example
/// assert_eq!(parse_motif_from_name("5xCAG").unwrap(), "CAG");
pub(crate) fn parse_motif_from_name(name: &str) -> Result<String> {
    let (repeat_str, motif) = name
        .split_once('x')
        .ok_or_else(|| Error::InvalidMsiBedMotif {
            motif: name.to_string(),
        })?;

    let repeat_count: usize = repeat_str.parse().map_err(|_| Error::InvalidMsiBedMotif {
        motif: name.to_string(),
    })?;

    // Validate motif contains only unambiguous DNA bases
    if !motif
        .chars()
        .all(|c| matches!(c, 'A' | 'C' | 'G' | 'T' | 'a' | 'c' | 'g' | 't'))
    {
        return Err(Error::InvalidMsiBedMotif {
            motif: name.to_string(),
        }
        .into());
    }

    match (repeat_count, motif.len()) {
        (0, _) => Err(Error::InvalidMsiBedMotif {
            motif: name.to_string(),
        }
        .into()),
        (_, 0) => Err(Error::InvalidMsiBedMotif {
            motif: name.to_string(), /* Note: If this is turned off in motif repeat check for msi module filtering should be turned on.*/
        }
        .into()),
        (_, _) => Ok(motif.to_string()),
    }
}

/*  NOTE: Other utility function(s) can be added that returns the repeat count
    of the motif togather with the motif or separately, if needed in the future,
    for extending the BedRegion struct.
*/

/// Parse a BED record into a BedRegion.
///
/// # Arguments
/// * `record` - BED record from bio::io::bed
///
/// # Returns
/// * `Ok(BedRegion)` - Parsed region with motif
/// * `Err` - Invalid coordinates or missing/invalid name
///
/// # Example
/// assert!(parse_bed_record(&record).is_ok());
pub(crate) fn parse_bed_record(record: &bed::Record) -> Result<BedRegion> {
    let chrom = record.chrom().to_string();
    let start = record.start();
    let end = record.end();

    if start >= end {
        return Err(Error::InvalidBedRecord {
            chrom: chrom.clone(),
            pos: start as i64,
            msg: "Invalid region coordinates: start >= end".to_string(),
        }
        .into());
    }

    let name = record.name().ok_or(Error::BedRecordMissingMotifName)?;
    let motif = parse_motif_from_name(name)?;

    Ok(BedRegion {
        chrom,
        start,
        end,
        motif,
    })
}

/// Validate BED file by checking first record.
///
/// Performs quick sanity check to ensure BED file is readable
/// and contains valid microsatellite region format.
///
/// # Arguments
/// * `bed_path` - Path to BED file
///
/// # Errors
/// Returns error if:
/// - File cannot be opened
/// - File is empty
/// - First record is invalid
///
/// # Example
/// assert!(validate_bed_file(&path).is_ok());
pub(crate) fn validate_bed_file(bed_path: &PathBuf) -> Result<()> {
    info!("Validating BED file: {}", bed_path.display());

    let mut bed_reader = bed::Reader::from_file(bed_path).context("Failed to open BED file")?;

    match bed_reader.records().next() {
        None => Err(Error::BedFileEmpty.into()),
        Some(Err(e)) => Err(e).context("Failed to read first BED record"),
        Some(Ok(record)) => {
            let first_region = parse_bed_record(&record)?;
            info!("BED file validated successfully");
            info!(
                "    First region: ({} {}-{} {}x{}) ms-status={}",
                first_region.chrom,
                first_region.start,
                first_region.end,
                first_region.motif_length(),
                first_region.motif,
                first_region.is_valid_motif(),
            );
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /* ============ BedRegion Tests ================== */

    #[test]
    fn test_bed_region_methods() {
        let region = BedRegion {
            chrom: "chr1".into(),
            start: 100,
            end: 200,
            motif: "CAG".into(),
        };

        assert_eq!(region.region_id(), "chr1:100-200");
        assert_eq!(region.motif_length(), 3);
        assert!(region.is_valid_motif());
    }

    #[test]
    fn test_bed_region_motif_validity() {
        // Valid: 1-6 bases
        for len in 1..=6 {
            let region = BedRegion {
                chrom: "chr1".into(),
                start: 0,
                end: 10,
                motif: "A".repeat(len),
            };
            assert!(region.is_valid_motif());
        }

        // Invalid: 0 or >6 bases
        let empty = BedRegion {
            chrom: "chr1".into(),
            start: 0,
            end: 10,
            motif: "".into(),
        };
        assert!(!empty.is_valid_motif());

        let too_long = BedRegion {
            chrom: "chr1".into(),
            start: 0,
            end: 10,
            motif: "AAAAAAA".into(),
        };
        assert!(!too_long.is_valid_motif());
    }

    /* ==== parse_motif_from_name tests ============== */

    #[test]
    fn test_parse_motif_valid() {
        assert_eq!(parse_motif_from_name("5xCAG").unwrap(), "CAG");
        assert_eq!(parse_motif_from_name("1xA").unwrap(), "A");
        assert_eq!(parse_motif_from_name("100xATCG").unwrap(), "ATCG");
    }

    #[test]
    fn test_parse_motif_invalid() {
        assert!(parse_motif_from_name("5CAG").is_err()); // Missing separator
        assert!(parse_motif_from_name("0xCAG").is_err()); // Zero count
        assert!(parse_motif_from_name("5x").is_err()); // Empty motif
        assert!(parse_motif_from_name("fivexCAG").is_err()); // Non-numeric count
    }

    #[test]
    fn test_parse_motif_invalid_chars() {
        // Motifs with ambiguous or invalid DNA letters should fail
        assert!(parse_motif_from_name("5xCAGN").is_err());
        assert!(parse_motif_from_name("3xXYZ").is_err());
        assert!(parse_motif_from_name("2xA-T").is_err());
        assert!(parse_motif_from_name("4x123").is_err());
    }

    /* ========= parse_bed_record tests ============== */

    /// Helper to create BED record from string
    fn bed_record_from_str(line: &str) -> bed::Record {
        let mut reader = bed::Reader::new(line.as_bytes());
        reader.records().next().unwrap().unwrap()
    }

    #[test]
    fn test_parse_bed_record_valid() {
        let record = bed_record_from_str("chr1\t100\t109\t3xCAG");

        let region = parse_bed_record(&record).unwrap();
        assert_eq!(region.chrom, "chr1");
        assert_eq!(region.start, 100);
        assert_eq!(region.end, 109);
        assert_eq!(region.motif, "CAG");

        /* TODO: Check This for MSI if we need to restrict small indels like 1xC */
        let record = bed_record_from_str("chr1\t100\t101\t1xC");
        assert!(parse_bed_record(&record).is_ok());
    }

    #[test]
    fn test_parse_bed_record_invalid_coords() {
        // start > end
        let record = bed_record_from_str("chr1\t109\t100\t3xCAG");
        assert!(parse_bed_record(&record).is_err());

        // start == end
        let record = bed_record_from_str("chr1\t100\t100\t1xC");
        assert!(parse_bed_record(&record).is_err());
    }

    #[test]
    fn test_parse_bed_record_missing_name() {
        let record = bed_record_from_str("chr1\t100\t200");
        assert!(parse_bed_record(&record).is_err());
    }

    /* ===== validate_bed_file tests ================= */

    #[test]
    fn test_validate_bed_file_valid() {
        let tmp = NamedTempFile::new().unwrap();
        writeln!(tmp.as_file(), "chr1\t100\t109\t3xCAG").unwrap();
        assert!(validate_bed_file(&tmp.path().to_path_buf()).is_ok());
    }

    #[test]
    fn test_validate_bed_file_empty() {
        let tmp = NamedTempFile::new().unwrap();
        assert!(validate_bed_file(&tmp.path().to_path_buf()).is_err());
    }

    #[test]
    fn test_validate_bed_file_invalid_record() {
        let tmp = NamedTempFile::new().unwrap();
        writeln!(tmp.as_file(), "chr1\t100\t200").unwrap(); // Missing name
        assert!(validate_bed_file(&tmp.path().to_path_buf()).is_err());
    }
}

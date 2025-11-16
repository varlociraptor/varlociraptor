//! ms_bed.rs
//!
//! Utilities for parsing and validating BED records used in MSI analysis.
//! The bed record format is expected to have the following fields:
//! chrom, start, end, region_id, (N)x(motif)
//! where (N)x(motif) represents a motif that is repeated N times.

use std::path::PathBuf;

use anyhow::{Context, Result};
use bio::io::bed;
use log::info;

use crate::errors::Error;

/// Represents a BED region entry with MSI motif information.
#[derive(Debug, Clone)]
pub(crate) struct BedRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub motif: String,
}

impl BedRegion {
    /// Returns formatted region ID like "chr1:100-200"
    pub(crate) fn region_id(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start, self.end)
    }

    /// Returns the motif length (derived from motif string).
    pub(crate) fn motif_length(&self) -> usize {
        self.motif.len()
    }

    /// Returns `true` if motif length is within [1,6].
    /// 1..=6 are valid motif lengths for MSI analysis.
    pub(crate) fn is_valid_motif(&self) -> bool {
        (1..=6).contains(&self.motif_length())
    }

    /*  NOTE: flexible is_valid_motif can be implemented if later need
        for MS regions to be greater than length 6.
    */
}

/// Parses a motif name string like "5xCAG" into motif and its length.
pub(crate) fn parse_motif_from_name(name: &str) -> Result<String> {
    let (repeat_str, motif) = name
        .split_once('x')
        .ok_or_else(|| Error::InvalidMsiBedMotif {
            motif: name.to_string(),
        })?;

    let repeat_count: usize = repeat_str.parse().map_err(|_| Error::InvalidMsiBedMotif {
        motif: name.to_string(),
    })?;

    match (repeat_count, motif.len()) {
        (0, _) => Err(Error::InvalidMsiBedMotif {
            motif: name.to_string(),
        }
        .into()),
        (_, 0) => Err(Error::InvalidMsiBedMotif {
            motif: name.to_string(), /* Note: If this is turned off in motif repeat check for msi mdoule filtering should be turned onn.*/
        }
        .into()),
        (_, _) => Ok(motif.to_string()),
    }
}

/*  NOTE: Other utility function(s) can be added that returns the repeat count
    of the motif togather with the motif or separately, if needed in the future,
    for extending the BedRegion struct.
*/

/// Parse a BED record into a BedRegion struct
pub(crate) fn parse_bed_record(record: &bed::Record) -> Result<BedRegion> {
    let chrom = record.chrom().to_string();
    let start = record.start();
    let end = record.end();

    if start >= end {
        return Err(Error::InvalidBCFRecord {
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

/// Quick sanity check(using single record) for BED file validity.
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

    #[test]
    fn test_parse_motif_ok() {
        let motif = parse_motif_from_name("5xCAG").unwrap();
        assert_eq!(motif, "CAG");
    }

    #[test]
    fn test_is_valid_motif_and_length() {
        let region = BedRegion {
            chrom: "chr1".into(),
            start: 0,
            end: 10,
            motif: "CAG".into(),
        };
        assert!(region.is_valid_motif());
        assert_eq!(region.motif_length(), 3);
    }
}

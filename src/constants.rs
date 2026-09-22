//! constants.rs
//!
//! This module defines global constants used across varlociraptor.
//! Centralizing constants here promotes maintainability and consistency across the codebase.
//!
//! Note: Constants should be defined under proper sections (e.g. MSI: CONSTANTS) with clear documentation on their purpose and usage.
//!
//! Sections of Constants included are:
//! 1. Generic: Constants - General constants used across the codebase (e.g. default values, thresholds).
//! 2. MSI: Constants - Microsatellite Instability related constants.
//! 3. Test: Constants - Floating-point comparison tolerances used only in test code.
//!

use rust_htslib::bcf::header::{TagLength, TagType};
use std::collections::HashSet;

/* =============== GENERIC: CONSTANTS ============== */

// Standard INFO fields always propagated in varlociraptor preprocessing.
// These fields are explicitly handled and must be omitted from aux_info.write()
// to prevent double-writing.
//
// Use standard_omit_aux_info!() macro instead of a constant for this string.

/// Floating-point tolerance for threshold/equality comparisons in production code
/// (e.g. checking whether a value is present in a set of expected values).
pub(crate) const EPSILON: f64 = 1e-9;

/// Minimum thread count for parallel tasks.
pub(crate) const MIN_THREAD_COUNT: usize = 1;

/* =============== MSI: CONSTANTS ================== */

/// Default MSI-High threshold (percentage).
///
/// Variants with MSI score >= this threshold are classified as MSI-High.
/// This is the standard clinical cutoff used in MSI analysis.
pub(crate) const MSI_DEFAULT_THRESHOLD: &str = "3.5";

/// Minimum allowed MSI threshold (percentage).
///
/// MSI thresholds must be positive values. Zero or negative thresholds
/// are not meaningful for MSI classification.
pub(crate) const MSI_MIN_THRESHOLD: f32 = 0.0;

/// Default AF thresholds for MSI evolution analysis.
/// This list is used to generate pseudo-time series of MSI evolution at
/// different allele frequency cutoffs. Also, this list is used to validate
/// the --distribution-af argument, as that needs to be from the list of
/// AF Thresholds.
pub(crate) const MSI_DEFAULT_AF_THRESHOLDS: &str = msi_default_af_thresholds!();

/// Sliding window size (bp) for regional MSI heatmap analysis.
pub(crate) const MSI_DEFAULT_SLIDING_WINDOW_SIZE: &str = "1000000";

/// INFO field tag for MSI region annotation.
/// Added to variants overlapping a microsatellite region.
pub const MSI_REGION_ID_TAG: &[u8] = b"REGION_ID";

/// INFO field tag for dummy deletion flag.
/// Set on synthetic records injected for MS regions without observed indels.
pub const MSI_DUMMY_TAG: &[u8] = b"MSI_DUMMY";

/// Default AF threshold for full distribution output.
/// Variants with AF >= this threshold are included in the full distribution output.
/// This AF threshold value is supported by the literature : https://pmc.ncbi.nlm.nih.gov/articles/PMC8172533/.
pub(crate) const MSI_DEFAULT_DISTRIBUTION_AF: &str = "0.05";

/// Default AF threshold for windowed heatmap analysis.
/// Variants with AF >= this threshold are included in the windowed heatmap analysis.
/// This AF threshold value is supported by the literature : https://pmc.ncbi.nlm.nih.gov/articles/PMC8172533/.
pub(crate) const MSI_DEFAULT_WINDOWED_AF: &str = "0.05";

/// Shape shared by varlociraptor-emitted per-ALT-allele Float
/// INFO/PROB_{event} fields.
pub(crate) const MSI_INFO_PROB_EVENT_FIELD_TYPE: TagType = TagType::Float;
pub(crate) const MSI_INFO_PROB_EVENT_FIELD_LENGTH: TagLength = TagLength::AltAlleles;

/// Shape shared by varlociraptor-emitted per-ALT-allele Float
/// FORMAT/AF field.
pub(crate) const MSI_FORMAT_AF_FIELD_TYPE: TagType = TagType::Float;
pub(crate) const MSI_FORMAT_AF_FIELD_LENGTH: TagLength = TagLength::AltAlleles;

lazy_static! {
    /// Full VCF header declaration for REGION_ID INFO field.
    /// Derived from MSI_REGION_ID_TAG.
    pub(crate) static ref MSI_REGION_ID_HEADER: String = format!(
        "##INFO=<ID={},Number=A,Type=String,Description=\"BED region ID for the overlapping microsatellite locus\">",
        std::str::from_utf8(MSI_REGION_ID_TAG).unwrap()
    );

    /// Full VCF header declaration for MSI_DUMMY INFO field.
    /// Derived from MSI_DUMMY_TAG.
    pub(crate) static ref MSI_DUMMY_HEADER: String = format!(
        "##INFO=<ID={},Number=0,Type=Flag,Description=\"Dummy deletion injected for MS region with no observed indel\">",
        std::str::from_utf8(MSI_DUMMY_TAG).unwrap()
    );

    /// INFO fields copied explicitly from input to output in MSI preprocessing.
    /// Derived from msi_omit_aux_info!() macro.
    pub(crate) static ref MSI_COPY_FIELDS: Vec<&'static str> =
        msi_omit_aux_info!()
            .split(", ")
            .collect();

    /// INFO fields omitted from aux_info.write() in MSI preprocessing.
    /// Combines standard propagated fields with MSI-specific output fields
    /// to prevent double-writing.
    pub(crate) static ref MSI_OMIT_AUX: HashSet<Vec<u8>> = {
        let msi_output_fields = [
            std::str::from_utf8(MSI_REGION_ID_TAG).unwrap(),
            std::str::from_utf8(MSI_DUMMY_TAG).unwrap(),
        ];
        msi_omit_aux_info!()
            .split(", ")
            .chain(msi_output_fields.iter().copied())
            .map(|s| s.as_bytes().to_vec())
            .collect()
    };
}

/* =============== TEST: CONSTANTS ================ */

#[cfg(test)]
pub(crate) mod test_constants {
    /// Tolerance for floating-point comparisons in tests - f64
    pub(crate) const TEST_EPSILON: f64 = 1e-6;

    /// Looser tolerance for floating-point comparisons (f64) where accumulated
    /// rounding error exceeds TEST_EPSILON - e.g. f32 storage round-trips,
    /// or multi-step probability conversions.
    pub(crate) const TEST_EPSILON_LOOSE: f64 = 1e-5;

    /// Tolerance for floating-point comparisons in tests - f32
    pub(crate) const TEST_EPSILON_F32: f32 = 1e-6;

    /// Looser tolerance for floating-point comparisons (f32) where accumulated
    /// rounding error exceeds TEST_EPSILON - e.g. f32 storage round-trips,
    /// or multi-step probability conversions.
    pub(crate) const TEST_EPSILON_LOOSE_F32: f32 = 1e-5;
}

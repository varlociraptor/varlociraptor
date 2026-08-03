//! mod.rs
//!
//! Microsatellite instability (MSI) estimation module.
//!
//! **Note: This feature is experimental.**
//!
//! This module provides:
//! 1. MSI configuration and validation
//! 2. Main pipeline orchestration
//!     - Coordination of intersection(intersection module),
//!     - Dynamic Programming analysis(dp_analysis module),
//!     - Output generation(output module)

mod dp_analysis;
mod intersection;
mod output;

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Result;
use log::info;

use crate::errors::Error;
use dp_analysis::OutputRequirements;

/* =============== CONSTANTS ====================== */

/// Default MSI-High threshold (percentage).
///
/// Variants with MSI score >= this threshold are classified as MSI-High.
/// This is the standard clinical cutoff used in MSI analysis.
pub const DEFAULT_MSI_THRESHOLD: &str = "3.5";

/// Minimum allowed MSI threshold (percentage).
///
/// MSI thresholds must be positive values. Zero or negative thresholds
/// are not meaningful for MSI classification.
pub const MIN_MSI_THRESHOLD: f64 = 0.0;

/* ================================================ */

/* ======== CLI CONFIGURATION ===================== */

/// Configuration for MSI estimation pipeline.
///
/// Contains all parameters needed for the MSI analysis workflow,
/// populated from CLI arguments(partially) and the file validation.
#[derive(Debug)]
pub(crate) struct MsiConfig {
    /// Path to BED file(should be sorted) with microsatellite regions.
    pub microsatellite_bed: PathBuf,
    /// Path to VCF/BCF file(should be sorted) with variant calls.
    pub calls: PathBuf,
    /// Number of threads (None = use rayon default).
    pub threads: Option<usize>,
    /// MSI-High classification threshold (percentage), default: 3.5.
    pub msi_threshold: f64,
    /// Sample names to process from VCF/BCF (after excluding user exclusions).
    pub samples: Vec<String>,
    /// Map of sample name to VCF header index (after excluding user exclusions).
    pub samples_index_map: HashMap<String, usize>,
    /// Whether the probabilities in the VCF are PHRED-scaled.
    pub is_phred: bool,
    /// Allele frequency thresholds to consider for AF evolution analysis
    /// when generating pseudotime outputs. (default: [1.0,0.8,0.6,0.4,0.2,0.0])
    /// If no pseudotime outputs are requested, this will be set to [0.0] to optimize computation.
    /// Note: This field is populated during CLI parsing and currently not validated for non [0-1] values
    /// as the this field is hidden constant set at CLI level. So future changes to expose this field
    /// to users should include validation for this field.
    pub af_thresholds: Vec<f64>,
    /// Output path for distribution plot (Vega-Lite JSON).
    pub plot_distribution: Option<PathBuf>,
    /// Output path for pseudotime plot (Vega-Lite JSON).
    pub plot_pseudotime: Option<PathBuf>,
    /// Output path for distribution data (TSV).
    pub data_distribution: Option<PathBuf>,
    /// Output path for pseudotime data (TSV).
    pub data_pseudotime: Option<PathBuf>,
}

impl MsiConfig {
    /// Validate the MSI configuration.
    /// Checks for valid MSI threshold and at least one output specified.
    pub fn validate(&self) -> Result<()> {
        if self.msi_threshold <= MIN_MSI_THRESHOLD {
            return Err(Error::MsiConfigThresholdInvalid {
                threshold: self.msi_threshold,
            }
            .into());
        }

        if self.plot_distribution.is_none()
            && self.plot_pseudotime.is_none()
            && self.data_distribution.is_none()
            && self.data_pseudotime.is_none()
        {
            return Err(Error::MsiConfigOutputMissing.into());
        }

        Ok(())
    }
}
/* ================================================ */

/* ============ MOD UTILITY ======================= */

/// Determine output requirements for computation optimization.
///
/// Analyzes which outputs are requested to determine:
/// - Whether to compute uncertainty bounds (needs_pseudotime)
/// - Whether to compute full probability distribution (needs_distribution)
///
/// This allows skipping expensive computations when not needed.
fn analyze_output_requirements(config: &MsiConfig) -> OutputRequirements {
    let needs_distribution =
        config.data_distribution.is_some() || config.plot_distribution.is_some();
    let needs_pseudotime = config.data_pseudotime.is_some() || config.plot_pseudotime.is_some();

    OutputRequirements {
        needs_pseudotime,
        needs_distribution,
    }
}

/* ================================================ */

/* ============ MAIN PIPELINE ===================== */

/// Main MSI estimation pipeline.
///
/// Orchestrates the complete MSI analysis workflow:
/// 1. **Intersection**: Stream BED regions against VCF variants
/// 2. **DP Analysis**: Compute probability distributions and MSI scores
/// 3. **Output**: Generate plots and data files
///
/// # Arguments
/// * `config` - Validated MSI configuration
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` on file I/O or processing errors
pub fn estimate_msi(config: MsiConfig) -> Result<()> {
    info!("----------------------------------------------");
    info!("Step 1: Config Stats");
    info!("----------------------------------------------");
    info!("BED file: {}", config.microsatellite_bed.display());
    info!("VCF/BCF file: {}", config.calls.display());
    info!("MSI threshold: {}", config.msi_threshold);
    info!("Samples included: {:?}", config.samples);

    info!("----------------------------------------------");
    info!("Step 2: Streaming Intersection");
    info!("----------------------------------------------");
    info!("Starting intersection of VCF/BCF and BED.");

    let (regions, total_regions) = intersection::intersect_streaming(
        &config.microsatellite_bed,
        &config.calls,
        &config.samples_index_map,
        config.is_phred,
    )?;

    if total_regions == 0 {
        info!("Terminating due to 0 valid regions.");
        return Ok(());
    }
    info!("Intersection finished successfully.");

    info!("----------------------------------------------");
    info!("Step 3: DP Analysis and Results Collection");
    info!("----------------------------------------------");

    // Determine computation requirements (single function)
    let output_requirements = analyze_output_requirements(&config);

    // Select AF thresholds based on output requirements for metrics optimization
    let af_thresholds: &[f64] = if output_requirements.needs_pseudotime {
        &config.af_thresholds
    } else {
        &[0.0]
    };

    let results = dp_analysis::run_af_evolution_analysis(
        &regions,
        total_regions,
        &config.samples,
        config.msi_threshold,
        af_thresholds,
        output_requirements,
        config.threads,
    )?;
    info!("Results generated successfully.");

    info!("----------------------------------------------");
    info!("Step 4(Final): Generating output(s)");
    info!("----------------------------------------------");

    if let Some(ref path) = &config.plot_distribution {
        output::generate_distribution_plot_spec(&results, path, config.msi_threshold)?;
        info!("Data Distribution Plot(Vega-Lite Json): {}", path.display());
    }

    if let Some(ref path) = config.plot_pseudotime {
        output::generate_pseudotime_plot_spec(&results, path, config.msi_threshold)?;
        info!("Pseudotime Plot(Vega-Lite Json): {}", path.display());
    }

    if let Some(ref path) = config.data_distribution {
        output::write_distribution_data(&results, path, config.msi_threshold)?;
        info!("Distribution Data(TSV): {}", path.display());
    }

    if let Some(ref path) = config.data_pseudotime {
        output::write_pseudotime_data(&results, path, config.msi_threshold)?;
        info!("Pseudotime Data(TSV): {}", path.display());
    }

    info!("==============================================");
    info!("MSI estimation complete");
    info!("==============================================");

    Ok(())
}

/* ================================================ */

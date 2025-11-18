//! mod.rs
//!
//! Microsatellite instability (MSI) estimation

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
pub const DEFAULT_MSI_THRESHOLD: &str = "3.5";
pub const MIN_MSI_THRESHOLD: f64 = 0.0;
/* ================================================ */

/* ======== CLI CONFIGURATION ===================== */
#[derive(Debug)]
pub(crate) struct MsiConfig {
    pub microsatellite_bed: PathBuf,
    pub calls: PathBuf,
    pub threads: Option<usize>,
    pub msi_threshold: f64,
    pub samples: Vec<String>,
    pub samples_index_map: HashMap<String, usize>,
    pub af_thresholds: Vec<f64>,
    pub plot_distribution: Option<PathBuf>,
    pub plot_pseudotime: Option<PathBuf>,
    pub data_distribution: Option<PathBuf>,
    pub data_pseudotime: Option<PathBuf>,
}

impl MsiConfig {
    pub fn validate(&self) -> Result<()> {
        if self.msi_threshold <= MIN_MSI_THRESHOLD {
            return Err(Error::InvalidMsiThreshold {
                threshold: self.msi_threshold,
            }
            .into());
        }

        if self.plot_distribution.is_none()
            && self.plot_pseudotime.is_none()
            && self.data_distribution.is_none()
            && self.data_pseudotime.is_none()
        {
            return Err(Error::NoMsiOutputSpecified.into());
        }

        Ok(())
    }
}
/* ================================================ */

/* ============ MOD UTILITY ======================= */
/// Determine what type of output requirements are required,
/// so that af thresholds can be setc accordingly and metrics
/// calculation can be optimized.
/// Returns output_requirements
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
/// Main Function to which CLI msi subcommand is dispatched to,
/// it runs the main pipeline orchestration.
///
/// This function coordinates the overall MSI estimation workflow by
/// streaming intersection, performing DP analysis, and generating outputs.
pub fn estimate_msi(config: MsiConfig) -> Result<()> {
    info!("**********************************************");
    info!("Config Stats: Streaming Intersection");
    info!("**********************************************");
    info!(" BED file: {}", config.microsatellite_bed.display());
    info!(" VCF/BCF file: {}", config.calls.display());
    info!(" MSI threshold: {}", config.msi_threshold);
    info!(" Samples included: {:?}", config.samples);

    info!("**********************************************");
    info!("Step 1: Streaming Intersection");
    info!("**********************************************");

    let (regions, total_regions) = intersection::intersect_streaming(
        &config.microsatellite_bed,
        &config.calls,
        &config.samples_index_map,
    )?;

    if total_regions == 0 {
        info!(" Terminating due to 0 valid regions.");
        return Ok(());
    }

    info!("**********************************************");
    info!("Step 2: DP Analysis and Results Collection");
    info!("**********************************************");

    // Determine computation requirements (single function)
    let output_requirements = analyze_output_requirements(&config);

    // Select AF thresholds based on output requirements for metrics optimization
    let af_thresholds: &[f64] = if output_requirements.needs_pseudotime {
        &config.af_thresholds
    } else {
        &[0.0]
    };
    info!(" AF Thresholds: {:?}", af_thresholds);

    let results = dp_analysis::run_af_evolution_analysis(
        &regions,
        total_regions,
        &config.samples,
        config.msi_threshold,
        af_thresholds,
        output_requirements,
        config.threads,
    )?;
    info!(" Results generated successfully.");

    info!("**********************************************");
    info!("Step 3(Final): Generating output(s)");
    info!("**********************************************");

    if let Some(ref path) = &config.plot_distribution {
        output::generate_distribution_plot_spec(&results, path, config.msi_threshold)?;
        info!(
            " Generated Data Distribution Plot(Vega-Lite Json): {}",
            path.display()
        );
    }

    if let Some(ref path) = config.plot_pseudotime {
        output::generate_pseudotime_plot_spec(&results, path, config.msi_threshold)?;
        info!(
            " Generated Pseudotime Plot(Vega-Lite Json): {}",
            path.display()
        );
    }

    if let Some(ref path) = config.data_distribution {
        output::write_distribution_data(&results, path, config.msi_threshold)?;
        info!(" Generated Distribution Data(TSV): {}", path.display());
    }

    if let Some(ref path) = config.data_pseudotime {
        output::write_pseudotime_data(&results, path, config.msi_threshold)?;
        info!(" Generated Pseudotime Data(TSV): {}", path.display());
    }

    info!("==============================================");
    info!("MSI estimation complete");
    info!("==============================================");

    Ok(())
}
/* ================================================ */

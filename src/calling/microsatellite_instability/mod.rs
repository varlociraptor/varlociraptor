//! mod.rs
//!
//! Microsatellite instability (MSI) calling module.
//! Calculates MSI scores from enriched, called VCF/BCF for a single sample.
//!
//! **Note: This feature is experimental.**
//!
//! This module provides:
//! Owns `MSIConfig` (CLI configuration) and `call_msi` (the entry point
//! that orchestrates extraction, analysis, and output). Declares these
//! child submodules:
//! 1. `extraction`  - extracts MSI regions and their variant/dummy status
//!    from the enriched, called VCF/BCF.
//! 2. `dp_analysis` - runs dynamic-programming AF-evolution analysis
//!    (global distribution/pseudotime, plus windowed heatmap analysis).
//! 3. `output`      - writes TSV data files and Vega-Lite plot specs for
//!    distribution, pseudotime, and heatmap outputs.
//! 4. `preprocessing` - contains the MSI preprocessing pipeline, which is
//!    used to generate the enriched VCF/BCF that is the input to this module.
//!

mod dp_analysis;
mod extraction;
mod output;
pub mod preprocessing;

use std::path::PathBuf;

use anyhow::{Context, Result};
use log::info;
use rust_htslib::bcf::{self, header::HeaderView, Read};

use crate::cli::{validate_thread_count, validate_vcf_file};
use crate::constants::{EPSILON, MSI_MIN_THRESHOLD};
use crate::errors::Error;
use crate::utils::bcf_utils::{
    get_sample_index, is_phred_scaled_from_path, validate_events_exist,
    validate_required_vcf_fields_msi, validate_samples_exist,
    validate_vcf_file as validate_ms_vcf_file,
};

use dp_analysis::AnalysisConfig;

/* ======== CLI CONFIGURATION ===================== */

/// Output requirements for optimizing calculations(conditional generation of metrics):
///
/// Determines which expensive computations are needed based on
/// which output files the user requested.
///     - needs_pseudotime: For pseudotime data
///     - needs_distribution: For distribution data
///     - needs_heatmap: For heatmap data
#[derive(Debug, Clone, Copy)]
struct OutputRequirements {
    /// Whether to compute uncertainty bounds (std dev, lower/upper)
    needs_pseudotime: bool,
    /// Whether to compute full probability distribution
    needs_distribution: bool,
    /// Whether to compute windowed heatmap analysis
    needs_heatmap: bool,
}

/// Configuration for MSI calling pipeline.
///
/// Contains all parameters needed for per-sample MSI analysis,
/// populated from CLI arguments.
#[derive(Debug)]
pub struct MSIConfig {
    /// Path to VCF/BCF file with variant calls.
    pub calls: PathBuf,
    /// Sample name to analyze from VCF/BCF.
    pub sample: String,
    /// Event names to combine for MSI probability (e.g., ["somatic_tumor", "high_vaf"].
    pub events: Vec<String>,
    /// MSI-High classification threshold (percentage), default: 3.5.
    pub msi_threshold: f32,
    /// Allele frequency thresholds to consider for AF evolution analysis
    /// when generating pseudotime outputs. (default: [1.0,0.8,0.6,0.4,0.2,0.0])
    /// If no pseudotime outputs are requested, this will be set to [0.0] to optimize computation.
    pub af_thresholds: Vec<f32>,
    /// Fixed AF at which the full P(K=k) distribution is reported.
    /// Must be present in `af_thresholds` - enforced in `validate()`.
    pub distribution_af: f32,
    /// Fixed AF used for windowed heatmap analysis. Independent of `af_thresholds`.
    pub windowed_af: f32,
    /// Sliding window size (bp) for regional MSI heatmap analysis (default: 1,000,000).
    /// Only used if --plot-heatmap or --data-heatmap specified.
    pub sliding_window: u64,
    /// Number of threads (None = use rayon default).
    pub threads: Option<usize>,
    /// Output path for distribution plot (Vega-Lite JSON).
    pub plot_distribution: Option<PathBuf>,
    /// Output path for pseudotime plot (Vega-Lite JSON).
    pub plot_pseudotime: Option<PathBuf>,
    /// Output path for heatmap plot (Vega-Lite JSON).
    pub plot_heatmap: Option<PathBuf>,
    /// Output path for distribution data (TSV).
    pub data_distribution: Option<PathBuf>,
    /// Output path for pseudotime data (TSV).
    pub data_pseudotime: Option<PathBuf>,
    /// Output path for heatmap data (TSV).
    pub data_heatmap: Option<PathBuf>,
    /// Whether the probabilities in the VCF are PHRED-scaled.
    pub is_phred: bool,
}

impl MSIConfig {
    /// Set default values for MSIConfig fields based on functional logic.
    /// Defaults Set:
    /// 1. Set is_phred based on events sepcification.
    pub fn set_defaults(&mut self) -> Result<()> {
        self.is_phred = is_phred_scaled_from_path(&self.calls)?;
        Ok(())
    }

    /// Validate the MSI configuration.
    ///
    /// Responsibilities:
    /// 1. `calls`: validates file extension, opens the file, checks it
    ///    contains at least one variant record, and that its header declares
    ///    FORMAT/AF with the correct type/shape (Float, one value per ALT allele).
    /// 2. `sample`: validates the calls file contains the specified sample.
    /// 3. `events`: validates the calls file contains the specified events in
    ///    the header, with each INFO/PROB_{EVENT} field having the correct
    ///    type/shape.
    /// 4. `msi_threshold`: checks for a valid MSI threshold.
    /// 5. `af_thresholds`: non-empty, and every entry in [0.0, 1.0]. Required
    ///    even though `--af-thresholds` is a `hidden` CLI flag - hidden only
    ///    suppresses it from `--help`, it does not prevent a caller (e.g. an
    ///    evaluation pipeline) from overriding it with an invalid value.
    /// 6. `distribution_af`: in [0.0, 1.0], and present in `af_thresholds`
    ///    - required because `calculate_msi_metrics` only populates `distribution
    ///    when `af_thresholds` contains distribution_af` exactly.
    /// 7. `windowed_af`: in [0.0, 1.0].
    /// 8. `sliding_window`: != 0 when a heatmap output is requested.
    /// 9. `threads`: if given, is >= `MIN_THREAD_COUNT`.
    /// 10. Output paths: checks at least one output is specified.
    pub fn validate(&self) -> Result<()> {
        // --- calls ---
        validate_vcf_file(&self.calls)?;
        let mut vcf =
            bcf::Reader::from_path(self.calls.as_path()).context("Failed to open VCF file")?;
        let header: HeaderView = vcf.header().clone();
        validate_ms_vcf_file(&mut vcf)?;
        validate_required_vcf_fields_msi(&header)?;

        // --- sample ---
        validate_samples_exist(&header, std::slice::from_ref(&self.sample))?;

        // --- events ---
        validate_events_exist(&header, &self.events)?;

        // --- msi_threshold ---
        if self.msi_threshold <= MSI_MIN_THRESHOLD {
            return Err(Error::MsiConfigThresholdInvalid {
                threshold: self.msi_threshold,
            }
            .into());
        }

        // --- af_thresholds ---
        if self.af_thresholds.is_empty() {
            return Err(Error::MsiConfigAfThresholdsEmpty.into());
        }

        for &af in &self.af_thresholds {
            if !(0.0..=1.0).contains(&af) {
                return Err(Error::MsiConfigAfThresholdsInvalid { threshold: af }.into());
            }
        }

        // --- distribution_af ---
        if !(0.0..=1.0).contains(&self.distribution_af) {
            return Err(Error::MsiConfigAfThresholdInvalid {
                field: "distribution_af",
                threshold: self.distribution_af,
            }
            .into());
        }

        let is_distribution_af = self
            .af_thresholds
            .iter()
            .any(|af| (af - self.distribution_af).abs() < EPSILON as f32);
        if !is_distribution_af {
            return Err(Error::MsiConfigDistributionAfMissing {
                distribution_af: self.distribution_af,
                af_thresholds: self.af_thresholds.clone(),
            }
            .into());
        }

        // --- windowed_af ---
        if !(0.0..=1.0).contains(&self.windowed_af) {
            return Err(Error::MsiConfigAfThresholdInvalid {
                field: "windowed_af",
                threshold: self.windowed_af,
            }
            .into());
        }

        // --- sliding_window ---
        let needs_heatmap = self.plot_heatmap.is_some() || self.data_heatmap.is_some();
        if needs_heatmap && self.sliding_window == 0 {
            return Err(Error::MsiConfigSlidingWindowInvalid {
                window_size: self.sliding_window,
            }
            .into());
        }

        // --- threads ---
        validate_thread_count(self.threads)?;

        // --- output paths ---
        if self.plot_distribution.is_none()
            && self.plot_pseudotime.is_none()
            && self.plot_heatmap.is_none()
            && self.data_distribution.is_none()
            && self.data_pseudotime.is_none()
            && self.data_heatmap.is_none()
        {
            return Err(Error::MsiConfigOutputMissing.into());
        }

        Ok(())
    }

    /// Log configuration details:
    ///   - Input file, sample, events
    ///   - Threshold and AF settings
    ///   - Requested outputs
    fn log_config(&self, output_req: &OutputRequirements) {
        info!("Input file: {}", self.calls.display());
        info!("Sample: {}", self.sample);
        info!("Events: {}", self.events.join(", "));
        info!("MSI threshold: {}", self.msi_threshold);
        if output_req.needs_pseudotime {
            info!(
                "AF thresholds: {}",
                self.af_thresholds
                    .iter()
                    .map(|af| af.to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
        if output_req.needs_distribution {
            info!("Distribution AF: {}", self.distribution_af);
        }
        if output_req.needs_heatmap {
            info!("Windowed AF: {}", self.windowed_af);
            info!("Sliding window: {} bp", self.sliding_window);
        }
        match self.threads {
            Some(t) => info!("Threads: {}", t),
            None => info!("Threads: default (rayon)"),
        }

        info!(
            "Probabilities are {} scaled",
            if self.is_phred { "PHRED" } else { "Linear" }
        );

        info!("Requested outputs:");
        if let Some(ref p) = self.data_distribution {
            info!("  - Distribution data (TSV): {}", p.display());
        }
        if let Some(ref p) = self.plot_distribution {
            info!("  - Distribution plot (JSON): {}", p.display());
        }
        if let Some(ref p) = self.data_pseudotime {
            info!("  - Pseudotime data (TSV): {}", p.display());
        }
        if let Some(ref p) = self.plot_pseudotime {
            info!("  - Pseudotime plot (JSON): {}", p.display());
        }
        if let Some(ref p) = self.data_heatmap {
            info!("  - Heatmap data (TSV): {}", p.display());
        }
        if let Some(ref p) = self.plot_heatmap {
            info!("  - Heatmap plot (JSON): {}", p.display());
        }
    }
}

/* ================================================ */

/// Orchestrates the MSI calling workflow: extract regions from a preprocessed +
/// called VCF/BCF, run AF-evolution analysis, then write requested outputs.
pub fn call_msi(config: MSIConfig) -> Result<()> {
    info!("----------------------------------------------");
    info!("Step 1: Configuration");
    info!("----------------------------------------------");

    let output_req = OutputRequirements {
        needs_pseudotime: config.plot_pseudotime.is_some() || config.data_pseudotime.is_some(),
        needs_distribution: config.plot_distribution.is_some()
            || config.data_distribution.is_some(),
        needs_heatmap: config.plot_heatmap.is_some() || config.data_heatmap.is_some(),
    };
    config.log_config(&output_req);

    info!("----------------------------------------------");
    info!("Step 2: Data Extraction");
    info!("----------------------------------------------");

    let (global_results, window_results) = {
        let mut vcf =
            bcf::Reader::from_path(&config.calls).context("Failed to open calls VCF/BCF")?;
        let header = vcf.header().clone();

        let sample_idx = get_sample_index(&header, &config.sample)?;

        let (regions, stats) = extraction::extract_regions(
            &mut vcf,
            &config.sample,
            sample_idx,
            &config.events,
            config.is_phred,
            output_req.needs_heatmap,
        )?;
        stats.log_stats();

        if stats.total_ms_regions == 0 {
            return Err(Error::MsiBedRegionsEmpty.into());
        }

        info!("----------------------------------------------");
        info!("Step 3: AF Evolution Analysis");
        info!("----------------------------------------------");

        let af_thresholds: Vec<f32> = if output_req.needs_pseudotime {
            config.af_thresholds.clone()
        } else if output_req.needs_distribution {
            vec![config.distribution_af]
        } else {
            vec![]
        };

        let analysis_config = AnalysisConfig {
            total_regions: stats.total_ms_regions,
            sample: &config.sample,
            msi_high_threshold: config.msi_threshold,
            af_thresholds,
            num_threads: config.threads,
            window_size: config.sliding_window,
            distribution_af: config.distribution_af,
            windowed_af: config.windowed_af,
        };

        dp_analysis::run_af_evolution_analysis(&regions, analysis_config, output_req)?
    };

    info!("----------------------------------------------");
    info!("Step 4: Output Generation");
    info!("----------------------------------------------");

    if let Some(ref path) = config.data_distribution {
        output::write_distribution_data(
            &global_results,
            &config.sample,
            path,
            config.msi_threshold,
            config.distribution_af,
        )?;
    }

    if let Some(ref path) = config.plot_distribution {
        output::generate_distribution_plot_spec(
            &global_results,
            &config.sample,
            path,
            config.msi_threshold,
            config.distribution_af,
        )?;
    }

    if let Some(ref path) = config.data_pseudotime {
        output::write_pseudotime_data(&global_results, &config.sample, path, config.msi_threshold)?;
    }

    if let Some(ref path) = config.plot_pseudotime {
        output::generate_pseudotime_plot_spec(
            &global_results,
            &config.sample,
            path,
            config.msi_threshold,
        )?;
    }

    if let Some(ref path) = config.data_heatmap {
        output::write_heatmap_data(
            &window_results,
            &config.sample,
            path,
            config.msi_threshold,
            config.windowed_af,
        )?;
    }

    if let Some(ref path) = config.plot_heatmap {
        output::generate_heatmap_plot_spec(
            &window_results,
            &config.sample,
            path,
            config.msi_threshold,
            config.windowed_af,
        )?;
    }

    info!("==============================================");
    info!("MSI calling complete");
    info!("==============================================");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::bcf_utils::tests::{create_test_vcf, TestVcfConfig};
    use std::path::Path;
    use tempfile::NamedTempFile;

    /// `create_test_vcf`'s temp file has no extension, but MSIConfig::validate()'s
    /// first check (crate::cli::validate_vcf_file) requires a .vcf-style filename.
    /// Copy the generated content into a properly-suffixed temp file.
    fn create_test_vcf_named(config: TestVcfConfig) -> (NamedTempFile, Vec<String>) {
        let (raw_tmp, sample_names) = create_test_vcf(config);
        let named_tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        std::fs::copy(raw_tmp.path(), named_tmp.path()).unwrap();
        (named_tmp, sample_names)
    }

    /// Build a config that passes every check, so each test can break exactly one field.
    fn valid_config(calls: &Path, sample: String) -> MSIConfig {
        MSIConfig {
            calls: calls.to_path_buf(),
            sample,
            events: vec!["somatic".to_string()], // matches PROB_SOMATIC in create_test_vcf's default output
            msi_threshold: 3.5,
            af_thresholds: vec![1.0, 0.5, 0.0],
            distribution_af: 0.5,
            windowed_af: 0.5,
            sliding_window: 1_000_000,
            threads: None,
            plot_distribution: Some(PathBuf::from("out.vl.json")),
            plot_pseudotime: None,
            plot_heatmap: None,
            data_distribution: None,
            data_pseudotime: None,
            data_heatmap: None,
            is_phred: false,
        }
    }

    #[test]
    fn test_validate_succeeds_with_valid_config() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let config = valid_config(tmp_vcf.path(), sample_names[0].clone());
        assert!(config.validate().is_ok());
    }

    #[test]
    fn test_validate_missing_sample() {
        let (tmp_vcf, _) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let config = valid_config(tmp_vcf.path(), "nonexistent_sample".to_string());

        let result = config.validate();
        assert!(result.is_err());
        assert!(format!("{}", result.unwrap_err()).contains("nonexistent_sample"));
    }

    #[test]
    fn test_validate_missing_event() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());
        config.events = vec!["not_a_real_event".to_string()];

        let result = config.validate();
        assert!(result.is_err());
        assert!(format!("{}", result.unwrap_err()).contains("not_a_real_event"));
    }

    #[test]
    fn test_validate_distribution_af_not_in_thresholds() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());
        config.af_thresholds = vec![1.0, 0.5, 0.0];
        config.distribution_af = 0.7; // not present in af_thresholds

        let result = config.validate();
        assert!(result.is_err());
        assert!(format!("{}", result.unwrap_err()).contains("distribution"));
    }

    #[test]
    fn test_validate_af_threshold_entry_out_of_range() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());
        config.af_thresholds = vec![1.5, 0.5]; // 1.5 is out of [0.0, 1.0]

        let result = config.validate();
        assert!(result.is_err());
        assert!(format!("{}", result.unwrap_err()).contains("1.5"));
    }

    #[test]
    fn test_validate_distribution_af_out_of_range() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());
        config.distribution_af = 1.5;

        let result = config.validate();
        assert!(result.is_err());
        assert!(format!("{}", result.unwrap_err()).contains("distribution_af"));
    }

    #[test]
    fn test_validate_msi_threshold_invalid() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());
        config.msi_threshold = 0.0; // must be > MIN_MSI_THRESHOLD (0.0)

        assert!(config.validate().is_err());
    }

    #[test]
    fn test_validate_af_thresholds_empty() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());
        config.af_thresholds = vec![];

        assert!(config.validate().is_err());
    }

    #[test]
    fn test_validate_windowed_af_out_of_range() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());
        config.windowed_af = 1.5; // out of [0.0, 1.0]

        let result = config.validate();
        assert!(result.is_err());
        assert!(format!("{}", result.unwrap_err()).contains("windowed_af"));
    }

    #[test]
    fn test_validate_sliding_window_zero_with_heatmap_requested() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());
        config.plot_heatmap = Some(PathBuf::from("heatmap.vl.json")); // triggers needs_heatmap
        config.sliding_window = 0;

        assert!(config.validate().is_err());
    }

    #[test]
    fn test_validate_threads_below_minimum() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());
        config.threads = Some(0); // MIN_THREAD_COUNT is 1

        assert!(config.validate().is_err());
    }

    #[test]
    fn test_validate_no_output_requested() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());
        config.plot_distribution = None; // valid_config's only output source, now cleared

        assert!(config.validate().is_err());
    }

    #[test]
    fn test_set_defaults_detects_linear_scale() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            use_phred: false,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());

        config.set_defaults().unwrap();
        assert!(!config.is_phred);
    }

    #[test]
    fn test_set_defaults_detects_phred_scale() {
        let (tmp_vcf, sample_names) = create_test_vcf_named(TestVcfConfig {
            num_samples: 1,
            use_phred: true,
            ..Default::default()
        });
        let mut config = valid_config(tmp_vcf.path(), sample_names[0].clone());

        config.set_defaults().unwrap();
        assert!(config.is_phred);
    }
}

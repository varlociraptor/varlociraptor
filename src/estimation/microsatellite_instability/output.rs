//! output.rs
//!
//! Output generation for MSI analysis.
//!
//! This module provides functions to generate four types of outputs:
//! 1. **Distribution data (TSV)**: Probability distribution P(K=k) at AF=0.0
//! 2. **Distribution plot (Vega-Lite JSON)**: Scatter plot of distribution
//! 3. **Pseudotime data (TSV)**: MSI score evolution across AF thresholds
//! 4. **Pseudotime plot (Vega-Lite JSON)**: Line plot with uncertainty bands
//!
//! All outputs are generated from `AfEvolutionResult` data computed by
//! the DP analysis module.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use serde_json::{json, Value};

use super::dp_analysis::AfEvolutionResult;

/* ============ Plotting Utils ==================== */

/// Create empty plot with informative message.
///
/// Uses embedded template from `templates/plots/msi_empty.json`.
/// Replaces "PLACEHOLDER_MESSAGE" with actual message.
///
/// # Arguments
/// * `path` - Output file path
/// * `message` - Text to display (e.g., "No distribution data available")
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation or JSON parsing fails
fn create_empty_plot(path: &Path, message: &str) -> Result<()> {
    let template = include_str!("../../../templates/plots/msi_empty.json");

    let mut spec: Value =
        serde_json::from_str(template).context("Failed to parse empty plot template")?;

    if let Value::Object(ref mut spec_obj) = spec {
        spec_obj["mark"]["text"] = json!(message);
    }

    let file = File::create(path).context("Failed to create empty plot file")?;
    let writer = BufWriter::new(file);

    serde_json::to_writer_pretty(writer, &spec).context("Failed to write empty plot")?;

    Ok(())
}

/// Load Vega-Lite template and inject data with threshold replacement.
///
/// # How It Works
/// 1. Parse template JSON from embedded string
/// 2. Inject data array into `data.values`
/// 3. Replace all threshold `datum` values in layers (X and Y encodings)
/// 4. Write pretty-printed JSON to file
///
/// # Threshold Replacement
/// Identifies threshold layers by:
/// - Mark type "rule" (threshold line)
/// - Mark type "text" with value "MSI Threshold" (label)
///
/// Updates both X and Y encodings where `datum` field exists.
///
/// # Arguments
/// * `data` - Data array as serde_json::Value
/// * `template` - Embedded template string
/// * `path` - Output file path
/// * `threshold` - MSI-High threshold (replaces all datum: 3.5)
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation, JSON parsing, or writing fails
fn write_plot(data: Value, template: &str, path: &Path, threshold: f64) -> Result<()> {
    let mut spec = serde_json::from_str(template).context("Failed to parse plot template")?;

    if let Value::Object(ref mut spec_obj) = spec {
        spec_obj["data"]["values"] = data;

        if let Some(layers) = spec_obj.get_mut("layer").and_then(|v| v.as_array_mut()) {
            for layer in layers {
                let mark_type = layer
                    .get("mark")
                    .and_then(|m| m.get("type"))
                    .and_then(|t| t.as_str());

                let is_threshold_layer = mark_type == Some("rule")
                    || (mark_type == Some("text")
                        && layer
                            .get("encoding")
                            .and_then(|e| e.get("text"))
                            .and_then(|t| t.get("value"))
                            .and_then(|v| v.as_str())
                            == Some("MSI Threshold"));

                if is_threshold_layer {
                    if let Some(x_encoding) = layer.get_mut("encoding").and_then(|e| e.get_mut("x"))
                    {
                        if x_encoding.get("datum").is_some() {
                            x_encoding["datum"] = json!(threshold);
                        }
                    }

                    if let Some(y_encoding) = layer.get_mut("encoding").and_then(|e| e.get_mut("y"))
                    {
                        if y_encoding.get("datum").is_some() {
                            y_encoding["datum"] = json!(threshold);
                        }
                    }
                }
            }
        }
    }

    let file = File::create(path).context("Failed to create plot file")?;
    let writer = BufWriter::new(file);

    serde_json::to_writer_pretty(writer, &spec).context("Failed to write plot JSON")?;

    Ok(())
}

/* ================================================ */

/* === Output Functions Type 1: Distribution ====== */

/// Write MSI probability distribution data to TSV file.
///
/// Outputs the complete probability distribution P(K=k) at AF=0.0 for each sample,
/// where k is the number of unstable microsatellite regions.
///
/// # Output Format
/// Tab-separated values with header:
/// ```text
/// sample  k   msi_score(threshold=3.5)    probability
/// tumor   0   0.00    0.420000
/// tumor   1   1.00    0.460000
/// tumor   2   2.00    0.120000
/// ```
///
/// # Data Availability
/// This function only writes data if:
/// - `output_req.needs_distribution = true` (set upstream in DP analysis)
/// - Results exist at AF=0.0 with distribution data
///
/// If no distribution data exists for any sample, writes header-only file.
///
/// # Arguments
/// * `results` - Nested HashMap: sample : AF threshold : MSI result
/// * `path` - Output TSV file path
/// * `msi_threshold` - MSI-High threshold for column header
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation or writing fails
///
/// # Example
/// ```no_run
/// write_distribution_data(&results, Path::new("dist.tsv"), 3.5)?;
/// ```
pub(super) fn write_distribution_data(
    results: &HashMap<String, HashMap<String, AfEvolutionResult>>,
    path: &Path,
    msi_threshold: f64,
) -> Result<()> {
    let file = File::create(path).context("Failed to create Distribution data TSV.")?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "sample\tk\tmsi_score(threshold={:.1})\tprobability",
        msi_threshold
    )?;

    for (sample, af_results) in results {
        let af_zero_result = match af_results.get("0") {
            Some(result) => result,
            None => continue,
        };

        let distribution = match &af_zero_result.distribution {
            Some(result) => result,
            None => continue,
        };

        for dp_result in distribution {
            writeln!(
                writer,
                "{}\t{}\t{:.2}\t{:.12}",
                sample, dp_result.k, dp_result.msi_score, dp_result.probability
            )?;
        }
    }

    writer.flush().context("Failed to flush distribution TSV")?;

    Ok(())
}

/// Generate Vega-Lite plot specification for MSI probability distribution.
///
/// Creates a scatter plot showing the probability distribution P(K=k) at AF=0.0,
/// with a vertical threshold line indicating the MSI-High cutoff.
///
/// Uses embedded template from `templates/plots/msi_distribution.json`.
///
/// # Plot Features
/// - **Scatter points**: MSI score (x-axis) vs posterior probability (y-axis)
/// - **Vertical threshold line**: Red dashed line at MSI-High cutoff
/// - **"MSI Threshold" label**: Text annotation at threshold position
/// - **Color by sample**: Automatic categorical colors for multiple samples
/// - **Interactive tooltips**: Sample name, MSI score, probability, k value
///
/// # Empty Case Handling
/// If no distribution data exists for any sample (no AF=0.0 results with
/// distribution), creates an empty plot with message "No distribution data available".
///
/// # Arguments
/// * `results` - Nested HashMap: sample : AF threshold : MSI result
/// * `path` - Output JSON file path (Vega-Lite specification)
/// * `msi_threshold` - MSI-High threshold for vertical line position
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation, template parsing, or writing fails
///
/// # Example
/// ```no_run
/// generate_distribution_plot_spec(&results, Path::new("dist.vl.json"), 3.5)?;
/// ```
pub(super) fn generate_distribution_plot_spec(
    results: &HashMap<String, HashMap<String, AfEvolutionResult>>,
    path: &Path,
    msi_threshold: f64,
) -> Result<()> {
    let mut data = Vec::new();

    for (sample, af_results) in results {
        let af_zero_result = match af_results.get("0") {
            Some(result) => result,
            None => continue,
        };

        let distribution = match &af_zero_result.distribution {
            Some(dist) => dist,
            None => continue,
        };

        for dp_result in distribution {
            data.push(json!({
                "sample": sample,
                "k": dp_result.k,
                "msi_score": dp_result.msi_score,
                "probability": dp_result.probability,
            }));
        }
    }

    if data.is_empty() {
        create_empty_plot(path, "No distribution data available")?;
        return Ok(());
    }

    write_plot(
        json!(data),
        include_str!("../../../templates/plots/msi_distribution.json"),
        path,
        msi_threshold,
    )?;

    Ok(())
}

/* ================================================ */

/* === Output Functions Type 2: Pseudotime ======== */

/// Write MSI pseudotime evolution trajectory data to TSV file.
///
/// Outputs MSI scores and uncertainty bounds across AF thresholds for each sample,
/// showing how MSI score evolves as variant detection threshold changes.
///
/// # Output Format
/// Tab-separated values with header:
/// **With uncertainty data (needs_pseudotime = true):**
/// sample    af_threshold    msi_score(threshold=3.5)    k_map    regions_with_variants    msi_status    uncertainty_lower    uncertainty_upper    map_std_dev
/// tumor     1.0             0.00                       0        20                       MSS          0.00                0.00                 0.000000
/// tumor     0.8             2.50                       2        35                       MSS          1.80                3.20                 0.700000
///
/// # Data Availability
/// Uncertainty columns (lower, upper, std_dev) are always included.
///
/// This function assumes uncertainty may exist because it's called when
/// `output_req.needs_pseudotime = true` (set upstream).
///
/// # Arguments
/// * `results` - Nested HashMap: sample : AF threshold : MSI result
/// * `path` - Output TSV file path
/// * `msi_threshold` - MSI-High threshold for column header
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation or writing fails
///
/// # Example
/// ```no_run
/// write_pseudotime_data(&results, Path::new("pseudo.tsv"), 3.5)?;
/// ```
pub(super) fn write_pseudotime_data(
    results: &HashMap<String, HashMap<String, AfEvolutionResult>>,
    path: &Path,
    msi_threshold: f64,
) -> Result<()> {
    let file = File::create(path).context("Failed to create Pseudotime data TSV.")?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "sample\taf_threshold\tmsi_score(threshold={:.1})\tk_map\tregions_with_variants\tmsi_status\tuncertainty_lower\tuncertainty_upper\tmap_std_dev",
        msi_threshold
    )?;

    for (sample, af_results) in results {
        let af_thresholds: Vec<String> = af_results.keys().cloned().collect();

        let mut af_pairs: Vec<(f64, String)> = af_thresholds
            .iter()
            .filter_map(|af_str| {
                af_str
                    .parse::<f64>()
                    .ok()
                    .map(|af_f64| (af_f64, af_str.clone()))
            })
            .collect();

        af_pairs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

        for (af_f64, af_str) in af_pairs {
            let result = match af_results.get(&af_str) {
                Some(r) => r,
                None => continue,
            };

            let lower = result
                .uncertainty_lower
                .map(|v| format!("{:.4}", v))
                .unwrap_or_else(|| "NA".to_string());

            let upper = result
                .uncertainty_upper
                .map(|v| format!("{:.4}", v))
                .unwrap_or_else(|| "NA".to_string());

            let std_dev = result
                .map_std_dev
                .map(|v| format!("{:.6}", v))
                .unwrap_or_else(|| "NA".to_string());

            writeln!(
                writer,
                "{}\t{:.1}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}",
                sample,
                af_f64,
                result.msi_score_map.unwrap(),
                result.k_map.unwrap(),
                result.regions_with_variants.unwrap(),
                result.msi_status.as_ref().unwrap(),
                lower,
                upper,
                std_dev
            )?;
        }
    }

    writer.flush().context("Failed to flush pseudotime TSV")?;
    Ok(())
}

/// Generate Vega-Lite plot specification for MSI pseudotime evolution.
///
/// Creates a line plot with uncertainty band showing MSI score trajectory
/// across AF thresholds, with a horizontal threshold line.
///
/// Uses embedded template from `templates/plots/msi_pseudotime.json`.
///
/// # Plot Features
/// - **Area band**: Shaded uncertainty range (lower to upper bound)
/// - **Line with points**: MSI score trajectory across AF thresholds
/// - **Horizontal threshold line**: Red dashed line at MSI-High cutoff
/// - **Threshold label**: Text annotation at threshold position
/// - **Color by sample**: Automatic categorical colors for multiple samples
/// - **Reversed X-axis**: 1.0 (left) to 0.0 (right), showing temporal flow
/// - **Interactive tooltips**: Sample name, AF threshold, MSI score, bounds
///
/// # Uncertainty Handling
/// If uncertainty bounds are missing for a data point, falls back to
/// using MSI score for both lower and upper bounds (zero-width band).
///
/// # Empty Case Handling
/// If no pseudotime data exists for any sample, creates an empty plot
/// with message "No pseudotime data available".
///
/// # Arguments
/// * `results` - Nested HashMap: sample : AF threshold : MSI result
/// * `path` - Output JSON file path (Vega-Lite specification)
/// * `msi_threshold` - MSI-High threshold for horizontal line position
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation, template parsing, or writing fails
///
/// # Example
/// ```ignore
/// generate_pseudotime_plot_spec(&results, Path::new("pseudo.vl.json"), 3.5)?;
/// ```
pub(super) fn generate_pseudotime_plot_spec(
    results: &HashMap<String, HashMap<String, AfEvolutionResult>>,
    path: &Path,
    msi_threshold: f64,
) -> Result<()> {
    let mut data = Vec::new();

    for (sample, af_results) in results {
        let mut af_pairs: Vec<(f64, String)> = af_results
            .keys()
            .filter_map(|af_str| {
                af_str
                    .parse::<f64>()
                    .ok()
                    .map(|af_f64| (af_f64, af_str.clone()))
            })
            .collect();

        af_pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        for (af_f64, af_str) in af_pairs {
            let result = match af_results.get(&af_str) {
                Some(r) => r,
                None => continue,
            };

            let msi_score = result.msi_score_map.unwrap();
            data.push(json!({
                "sample": sample,
                "af_threshold": af_f64,
                "msi_score": result.msi_score_map,
                "lower_bound": result.uncertainty_lower.unwrap_or(msi_score),
                "upper_bound": result.uncertainty_upper.unwrap_or(msi_score),
            }));
        }
    }

    if data.is_empty() {
        create_empty_plot(path, "No pseudotime data available")?;
        return Ok(());
    }

    write_plot(
        json!(data),
        include_str!("../../../templates/plots/msi_pseudotime.json"),
        path,
        msi_threshold,
    )
}

/* ================================================ */

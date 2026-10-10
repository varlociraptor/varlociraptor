//! output.rs
//!
//! Output generation for MSI analysis.
//!
//! This module provides functions to generate six types of outputs:
//! 1. **Distribution data (TSV)**: Probability distribution P(K=k) at AF=0.0
//! 2. **Distribution plot (Vega-Lite JSON)**: Scatter plot of distribution
//! 3. **Pseudotime data (TSV)**: MSI score evolution across AF thresholds
//! 4. **Pseudotime plot (Vega-Lite JSON)**: Line plot with uncertainty bands
//! 5. **Heatmap data (TSV)**: Windowed MSI scores across the genome
//! 6. **Heatmap plot (Vega-Lite JSON)**: Spatial rect heatmap by chromosome
//!
//! Types 1-4 use `AfEvolutionResult` from DP analysis.
//! Types 5-6 use `WindowResult` from windowed analysis.
//!
//! # Output Requirements
//! Each function is called only when the corresponding `OutputRequirements`
//! flag is set, checked in upstream.
//!
//! NOTE: An empty-plot placeholder function (create_empty_plot) was removed
//! as unreachable dead code - available to copy back if ever needed:
//! https://github.com/rohan-ibn-tariq/varlociraptor/blob/56278ba1f36f8a89046c3bc9481d502ab7e0b377/src/estimation/microsatellite_instability/output.rs#L26
//!

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use serde_json::{json, Value};

use super::dp_analysis::{AfEvolutionResult, WindowResult};
use crate::utils::genomics::classify_msi_status;

/* ============ Data Structures =================== */

/// Determines which template mutations `write_plot` applies
/// when injecting threshold values into a Vega-Lite specification.
enum PlotType {
    /// Updates datum in rule/text layers - works for both vertical (distribution)
    /// and horizontal (pseudotime) threshold lines since the layer detection
    /// finds whichever axis has a datum field.
    WithThresholdLine,
    /// Updates color scale domain midpoint for diverging coloring.
    Heatmap,
}

/* ================================================ */

/* ============ Plotting Utils ==================== */

/// Load Vega-Lite template, inject data, apply threshold mutation,
/// and write pretty-printed JSON to file.
///
/// # How It Works
/// 1. Parse template JSON from embedded string
/// 2. Inject data array into `data.values`
/// 3. Apply plot-type-specific threshold mutation (see below)
/// 4. Write pretty-printed JSON to file
///
/// # Plot Type Behavior
/// - `WithThresholdLine`: Finds rule and text layers by mark type,
///   then updates whichever of X or Y encoding has a `datum` field.
///   Distribution uses X datum (vertical line), Pseudotime uses Y datum
///   (horizontal line) - detection handles both without branching.
/// - `Heatmap`: Updates `encoding.color.scale.domain[1]` to set the
///   diverging color scale midpoint at the MSI-High threshold.
///
/// # Arguments
/// * `data`      - Data array to inject into `data.values`
/// * `template`  - Embedded template string
/// * `path`      - Output file path
/// * `threshold` - MSI-High threshold value
/// * `plot_type` - Controls which template mutation is applied
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation, JSON parsing, or writing fails
fn write_plot(
    data: Value,
    template: &str,
    path: &Path,
    threshold: f32,
    plot_type: PlotType,
    title: &str,
) -> Result<()> {
    let mut spec: Value =
        serde_json::from_str(template).context("Failed to parse plot template")?;

    if let Value::Object(ref mut spec_obj) = spec {
        spec_obj["data"]["values"] = data;
        spec_obj["title"] = json!(title);

        // This substitution depends on hardcoded template details (a "rule" mark,
        // "MSI Threshold" text, a 3-value color list) - changing the templates
        // can silently break the plotting without this code (or output) noticing.
        match plot_type {
            PlotType::WithThresholdLine => {
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

                        if !is_threshold_layer {
                            continue;
                        }

                        // Update whichever axis has a datum field
                        // Distribution uses X, Pseudotime uses Y -
                        // detection handles both without branching
                        for axis in ["x", "y"] {
                            if let Some(enc) =
                                layer.get_mut("encoding").and_then(|e| e.get_mut(axis))
                            {
                                if enc.get("datum").is_some() {
                                    enc["datum"] = json!(threshold);
                                }
                            }
                        }
                    }
                }
            }

            PlotType::Heatmap => {
                if let Some(domain) = spec_obj
                    .get_mut("encoding")
                    .and_then(|e| e.get_mut("color"))
                    .and_then(|c| c.get_mut("scale"))
                    .and_then(|s| s.get_mut("domain"))
                    .and_then(|d| d.as_array_mut())
                {
                    if domain.len() == 3 {
                        domain[1] = json!(threshold);
                    }
                }
            }
        }
    }

    let file = File::create(path).context("Failed to create plot file")?;
    let mut writer = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, &spec).context("Failed to write plot JSON")?;
    writer.flush().context("Failed to flush plot file")?;

    Ok(())
}

/* ================================================ */

/* === Output Functions Type 1: Distribution ====== */

/// Write MSI probability distribution data to TSV file.
///
/// Outputs the complete probability distribution P(K=k) at AF=0.0,
/// where k is the number of unstable microsatellite regions.
///
/// # Output Format
/// Tab-separated values with header:
/// ```text
/// sample  distribution_af  k   msi_score(threshold=3.5)    probability
/// tumor   0.00             0   0.00    0.420000
/// tumor   0.00             1   1.00    0.460000
/// tumor   0.00             2   2.00    0.120000
/// ```
///
/// # Arguments
/// * `results`       - Flat HashMap keyed by AF threshold string (e.g. "0.00")
/// * `sample`        - Sample name written to the first column of each row
/// * `path`          - Output TSV file path
/// * `msi_threshold` - MSI-High threshold for column header label
/// * `distribution_af` - Allele frequency threshold for which distribution is written
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation or writing fails
///
/// # Example
/// ```ignore
/// write_distribution_data(&results, "sample", Path::new("dist.tsv"), 3.5, 0.0)?;
/// ```
pub(super) fn write_distribution_data(
    results: &HashMap<String, AfEvolutionResult>,
    sample: &str,
    path: &Path,
    msi_threshold: f32,
    distribution_af: f32,
) -> Result<()> {
    let file = File::create(path).context("Failed to create Distribution data TSV.")?;
    let mut writer = BufWriter::new(file);
    let af_key = format!("{:.2}", distribution_af);

    writeln!(
        writer,
        "sample\tdistribution_af\tk\tmsi_score(threshold={:.1})\tprobability",
        msi_threshold
    )?;

    if let Some(distribution) = results.get(&af_key).and_then(|r| r.distribution.as_ref()) {
        for dp_result in distribution {
            writeln!(
                writer,
                "{}\t{:.2}\t{}\t{:.2}\t{:.12}",
                sample, distribution_af, dp_result.k, dp_result.msi_score, dp_result.probability
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
/// - **Interactive tooltips**: MSI score, probability, k value
///
/// # Arguments
/// * `results`       - Flat HashMap keyed by AF threshold string (e.g. "0.00")
/// * `sample`        - Sample name for plot title.
/// * `path`          - Output JSON file path (Vega-Lite specification)
/// * `msi_threshold` - MSI-High threshold for vertical line position
/// * `distribution_af` - AF threshold this distribution was computed at
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation, template parsing, or writing fails
///
/// # Example
/// ```ignore
/// generate_distribution_plot_spec(&results, "tumor", Path::new("dist.vl.json"), 3.5, 0.0)?;
/// ```
pub(super) fn generate_distribution_plot_spec(
    results: &HashMap<String, AfEvolutionResult>,
    sample: &str,
    path: &Path,
    msi_threshold: f32,
    distribution_af: f32,
) -> Result<()> {
    let mut data = Vec::new();
    let af_key = format!("{:.2}", distribution_af);

    if let Some(distribution) = results.get(&af_key).and_then(|r| r.distribution.as_ref()) {
        for dp_result in distribution {
            data.push(json!({
                "k": dp_result.k,
                "msi_score": dp_result.msi_score,
                "probability": dp_result.probability,
            }));
        }
    }

    write_plot(
        json!(data),
        include_str!("../../../templates/plots/msi_distribution.json"),
        path,
        msi_threshold,
        PlotType::WithThresholdLine,
        &format!(
            "MSI Score Distribution — Sample: {} (AF ≥ {:.2})",
            sample, distribution_af
        ),
    )?;

    Ok(())
}

/* ================================================ */

/* === Output Functions Type 2: Pseudotime ======== */

/// Write MSI pseudotime evolution trajectory data to TSV file.
///
/// Outputs MSI scores and uncertainty bounds across AF thresholds,
/// showing how MSI score evolves as variant detection threshold changes.
///
/// # Output Format
/// Tab-separated values, sorted highest to lowest AF threshold.
/// ```text
/// sample    af_threshold    total_regions    msi_score(threshold=3.5)    k_map    regions_af_pass    msi_status    uncertainty_lower    uncertainty_upper    map_std_dev
/// tumor     1.0             40               0.00                       0        20                 MSS          0.00                0.00                 0.000000
/// tumor     0.8             40               2.50                       2        35                 MSS          1.80                3.20                 0.700000
/// ```
///
/// # Arguments
/// * `results`       - Flat HashMap keyed by AF threshold string
/// * `sample`        - Sample name written to the first column of each row
/// * `path`          - Output TSV file path
/// * `msi_threshold` - MSI-High threshold for column header
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation or writing fails
///
/// # Example
/// ```ignore
/// write_pseudotime_data(&results, "tumor", Path::new("pseudo.tsv"), 3.5)?;
/// ```
pub(super) fn write_pseudotime_data(
    results: &HashMap<String, AfEvolutionResult>,
    sample: &str,
    path: &Path,
    msi_threshold: f32,
) -> Result<()> {
    let file = File::create(path).context("Failed to create Pseudotime data TSV.")?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "sample\taf_threshold\ttotal_regions\tmsi_score(threshold={:.1})\tk_map\tregions_af_pass\tmsi_status\tuncertainty_lower\tuncertainty_upper\tmap_std_dev",
        msi_threshold
    )?;

    let mut af_pairs: Vec<(f32, String)> = results
        .keys()
        .filter_map(|af_str| {
            af_str
                .parse::<f32>()
                .ok()
                .map(|af_f32| (af_f32, af_str.clone()))
        })
        .collect();

    af_pairs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    for (af_f32, af_str) in af_pairs {
        let result = match results.get(&af_str) {
            Some(r) => r,
            None => continue,
        };

        let msi_score_str = result
            .msi_score_map
            .map(|v| format!("{:.2}", v))
            .unwrap_or_else(|| "NA".to_string());
        let k_map_str = result
            .k_map
            .map(|v| v.to_string())
            .unwrap_or_else(|| "NA".to_string());
        let regions_str = result
            .regions_af_pass
            .map(|v| v.to_string())
            .unwrap_or_else(|| "NA".to_string());
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
            "{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            sample,
            af_f32,
            result.total_regions,
            msi_score_str,
            k_map_str,
            regions_str,
            result
                .msi_score_map
                .map(|s| classify_msi_status(s, msi_threshold).as_str())
                .unwrap_or("NA"),
            lower,
            upper,
            std_dev
        )?;
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
/// - **Reversed X-axis**: 1.0 (left) to 0.0 (right), showing temporal flow
/// - **Interactive tooltips**: AF threshold, MSI score, bounds
///
/// # Uncertainty Handling
/// If uncertainty bounds are absent, MSI score is used for both bounds.
///
/// # Arguments
/// * `results`       - Flat HashMap keyed by AF threshold string
/// * `sample`        - Sample name for plot title.
/// * `path`          - Output JSON file path (Vega-Lite specification)
/// * `msi_threshold` - MSI-High threshold for horizontal line position
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation, template parsing, or writing fails
///
/// # Example
/// ```ignore
/// generate_pseudotime_plot_spec(&results, "tumor", Path::new("pseudo.vl.json"), 3.5)?;
/// ```
pub(super) fn generate_pseudotime_plot_spec(
    results: &HashMap<String, AfEvolutionResult>,
    sample: &str,
    path: &Path,
    msi_threshold: f32,
) -> Result<()> {
    let mut data = Vec::new();

    let mut af_pairs: Vec<(f32, String)> = results
        .keys()
        .filter_map(|af_str| {
            af_str
                .parse::<f32>()
                .ok()
                .map(|af_f32| (af_f32, af_str.clone()))
        })
        .collect();

    af_pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    for (af_f32, af_str) in &af_pairs {
        let result = match results.get(af_str) {
            Some(r) => r,
            None => continue,
        };

        let msi_score = match result.msi_score_map {
            Some(v) => v,
            None => continue, // no computed score for this AF threshold - skip the point
        };

        data.push(json!({
            "af_threshold": af_f32,
            "msi_score": msi_score,
            "lower_bound": result.uncertainty_lower.unwrap_or(msi_score),
            "upper_bound": result.uncertainty_upper.unwrap_or(msi_score),
        }));
    }

    write_plot(
        json!(data),
        include_str!("../../../templates/plots/msi_pseudotime.json"),
        path,
        msi_threshold,
        PlotType::WithThresholdLine,
        &format!("MSI Evolution Across AF Thresholds — Sample: {}", sample),
    )
}

/* ================================================ */

/* === Output Functions Type 3: Heatmap =========== */

/// Write windowed MSI heatmap data to TSV file.
///
/// Each row represents one genomic window with its MSI score,
/// posterior probability, and MSI-High classification.
///
/// # Output Format
/// ```text
/// sample  windowed_af  chrom  window_start  window_end  msi_score(threshold=3.5)  posterior_probability  regions_in_window  msi_status
/// tumor   0.05         chr1   0             1000000     2.5000                    0.820000               40                 MSS
/// ```
///
/// # Arguments
/// * `results`       - Windowed analysis results from `run_windowed_analysis`
/// * `sample`        - Sample name for the first column
/// * `path`          - Output TSV file path
/// * `msi_threshold` - MSI-High threshold for column header and classification
/// * `windowed_af`   - Fixed AF threshold windowed analysis was computed at
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation or writing fails
///
/// # Example
/// ```ignore
/// write_heatmap_data(&results, "tumor", Path::new("heatmap.tsv"), 3.5, 0.05)?;
/// ```
pub(super) fn write_heatmap_data(
    results: &[WindowResult],
    sample: &str,
    path: &Path,
    msi_threshold: f32,
    windowed_af: f32,
) -> Result<()> {
    let file = File::create(path).context("Failed to create heatmap TSV")?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "sample\twindowed_af\tchrom\twindow_start\twindow_end\tmsi_score(threshold={:.1})\tposterior_probability\tregions_in_window\tmsi_status",
        msi_threshold
    )?;

    for r in results {
        writeln!(
            writer,
            "{}\t{:.2}\t{}\t{}\t{}\t{:.4}\t{:.6}\t{}\t{}",
            sample,
            windowed_af,
            r.chrom,
            r.window_start,
            r.window_end,
            r.msi_score,
            r.posterior_probability,
            r.regions_in_window,
            classify_msi_status(r.msi_score, msi_threshold),
        )?;
    }

    writer.flush().context("Failed to flush heatmap TSV")?;
    Ok(())
}

/// Generate Vega-Lite plot specification for windowed MSI heatmap.
///
/// Creates a diverging-color rect heatmap where:
/// - X-axis: genomic window position (start to end)
/// - Y-axis: chromosome (ordinal, genomic order preserved via sort:null)
/// - Color:  MSI score - blue (stable) through soft yellow (threshold) to red (unstable)
/// - Opacity: posterior probability (confidence in the score)
/// - Tooltip: includes MSI-High/MSS classification at the given threshold
///
/// Uses `write_plot` with embedded template `templates/plots/msi_heatmap.json`.
/// `write_plot` dynamically updates the color scale midpoint to match
/// `msi_threshold` - the template uses a 3-point domain [0, threshold, 100].
///
/// # Arguments
/// * `results`       - Windowed analysis results from `run_windowed_analysis`
/// * `sample`        - Sample name for title
/// * `path`          - Output JSON file path
/// * `msi_threshold` - MSI-High threshold for color scale midpoint and classification
/// * `windowed_af`   - Fixed AF threshold windowed analysis was computed at
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if file creation, template parsing, or writing fails
///
/// # Example
/// ```ignore
/// generate_heatmap_plot_spec(&results, "tumor", Path::new("heatmap.vl.json"), 3.5, 0.05)?;
/// ```
pub(super) fn generate_heatmap_plot_spec(
    results: &[WindowResult],
    sample: &str,
    path: &Path,
    msi_threshold: f32,
    windowed_af: f32,
) -> Result<()> {
    let data: Vec<Value> = results
        .iter()
        .map(|r| {
            json!({
                "chrom": r.chrom,
                "window_start": r.window_start,
                "window_end": r.window_end,
                "msi_score": r.msi_score,
                "posterior_probability": r.posterior_probability,
                "regions_in_window": r.regions_in_window,
                "msi_status": classify_msi_status(r.msi_score, msi_threshold).as_str(),
            })
        })
        .collect();

    write_plot(
        json!(data),
        include_str!("../../../templates/plots/msi_heatmap.json"),
        path,
        msi_threshold,
        PlotType::Heatmap,
        &format!(
            "MSI Spatial Heatmap — Sample: {} (AF ≥ {:.2})",
            sample, windowed_af
        ),
    )
}

/* ================================================ */

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::NamedTempFile;

    use crate::calling::microsatellite_instability::dp_analysis::{
        AfEvolutionResult, DpResult, WindowResult,
    };
    use crate::utils::genomics::MsiStatus;

    /* ====== Test helpers ======================== */

    fn make_distribution_result() -> HashMap<String, AfEvolutionResult> {
        let mut m = HashMap::new();
        m.insert(
            "0.00".to_string(),
            AfEvolutionResult {
                sample: "tumor".to_string(),
                af_threshold: 0.0,
                total_regions: 40,
                k_map: None,
                msi_score_map: None,
                regions_af_pass: None,
                uncertainty_lower: None,
                uncertainty_upper: None,
                map_std_dev: None,
                distribution: Some(vec![
                    DpResult {
                        k: 0,
                        msi_score: 0.0,
                        probability: 0.6,
                    },
                    DpResult {
                        k: 1,
                        msi_score: 50.0,
                        probability: 0.4,
                    },
                ]),
            },
        );
        m
    }

    fn make_pseudotime_result() -> HashMap<String, AfEvolutionResult> {
        let mut m = HashMap::new();
        m.insert(
            "0.00".to_string(),
            AfEvolutionResult {
                sample: "tumor".to_string(),
                af_threshold: 0.0,
                total_regions: 40,
                k_map: Some(2),
                msi_score_map: Some(2.0),
                regions_af_pass: Some(10),
                uncertainty_lower: Some(1.0),
                uncertainty_upper: Some(3.0),
                map_std_dev: Some(0.5),
                distribution: None,
            },
        );
        m
    }

    fn make_window_results() -> Vec<WindowResult> {
        vec![WindowResult {
            chrom: "chr1".to_string(),
            window_start: 0,
            window_end: 1_000_000,
            msi_score: 2.5,
            posterior_probability: 0.8,
            regions_in_window: 40,
        }]
    }

    /* ====== Distribution TSV =================== */

    #[test]
    fn test_write_distribution_data_row_written() {
        let tmp = NamedTempFile::new().unwrap();
        write_distribution_data(&make_distribution_result(), "tumor", tmp.path(), 3.5, 0.0)
            .unwrap();

        let content = fs::read_to_string(tmp.path()).unwrap();
        assert_eq!(content.lines().count(), 3); // header + 2 rows

        let first_row = content.lines().nth(1).unwrap();
        assert!(first_row.starts_with("tumor"));
        assert!(first_row.contains("0.00")); // msi_score k=0
        assert!(first_row.contains("0.6")); // probability k=0
    }

    #[test]
    fn test_generate_distribution_plot_threshold_injected() {
        let tmp = NamedTempFile::new().unwrap();
        generate_distribution_plot_spec(&make_distribution_result(), "tumor", tmp.path(), 7.5, 0.0)
            .unwrap();

        let content = fs::read_to_string(tmp.path()).unwrap();
        let spec: Value = serde_json::from_str(&content).unwrap();

        let layers = spec["layer"].as_array().expect("plot should have layers");

        let threshold_layers: Vec<&Value> = layers
            .iter()
            .filter(|layer| {
                let mark_type = layer["mark"]["type"].as_str();
                mark_type == Some("rule")
                    || layer["encoding"]["text"]["value"].as_str() == Some("MSI Threshold")
            })
            .collect();

        assert_eq!(
            threshold_layers.len(),
            2,
            "expected both the rule layer and the text layer to match"
        );

        for layer in &threshold_layers {
            assert_eq!(
                layer["encoding"]["x"]["datum"].as_f64(),
                Some(7.5),
                "every threshold layer's x.datum should be updated"
            );
        }
    }

    /* ====== Pseudotime TSV ===================== */

    #[test]
    fn test_write_pseudotime_data_row_written() {
        let tmp = NamedTempFile::new().unwrap();
        write_pseudotime_data(&make_pseudotime_result(), "tumor", tmp.path(), 3.5).unwrap();

        let content = fs::read_to_string(tmp.path()).unwrap();
        assert_eq!(content.lines().count(), 2); // header + one row

        let row = content.lines().nth(1).unwrap();
        assert!(row.starts_with("tumor"));
        assert!(row.contains("0.00")); // af_threshold
        assert!(row.contains("2.00")); // msi_score
        assert!(row.contains(MsiStatus::Stable.as_str())); // MSI status
    }

    #[test]
    fn test_generate_pseudotime_plot_threshold_injected() {
        let tmp = NamedTempFile::new().unwrap();
        generate_pseudotime_plot_spec(&make_pseudotime_result(), "tumor", tmp.path(), 7.5).unwrap();

        let content = fs::read_to_string(tmp.path()).unwrap();
        let spec: Value = serde_json::from_str(&content).unwrap();

        let layers = spec["layer"].as_array().expect("plot should have layers");

        let threshold_layers: Vec<&Value> = layers
            .iter()
            .filter(|layer| {
                let mark_type = layer["mark"]["type"].as_str();
                mark_type == Some("rule")
                    || layer["encoding"]["text"]["value"].as_str() == Some("MSI Threshold")
            })
            .collect();

        assert_eq!(
            threshold_layers.len(),
            2,
            "expected both the rule layer and the text layer to match"
        );

        for layer in &threshold_layers {
            assert_eq!(
                layer["encoding"]["y"]["datum"].as_f64(),
                Some(7.5),
                "every threshold layer's y.datum should be updated"
            );
        }
    }

    /* ====== Heatmap TSV ======================== */

    #[test]
    fn test_write_heatmap_data_row_written() {
        let tmp = NamedTempFile::new().unwrap();
        write_heatmap_data(&make_window_results(), "tumor", tmp.path(), 3.5, 0.05).unwrap();

        let content = fs::read_to_string(tmp.path()).unwrap();
        assert_eq!(content.lines().count(), 2); // header + one row

        let row = content.lines().nth(1).unwrap();
        assert!(row.starts_with("tumor"));
        assert!(row.contains("chr1"));
        assert!(row.contains(MsiStatus::Stable.as_str())); // MSI status
        assert!(row.contains("0.8000")); // posterior
    }

    #[test]
    fn test_generate_heatmap_plot_threshold_injected() {
        let tmp = NamedTempFile::new().unwrap();
        generate_heatmap_plot_spec(&make_window_results(), "tumor", tmp.path(), 7.5, 0.05).unwrap();

        let content = fs::read_to_string(tmp.path()).unwrap();
        let spec: Value = serde_json::from_str(&content).unwrap();

        let domain = spec["encoding"]["color"]["scale"]["domain"]
            .as_array()
            .expect("heatmap should have a color domain");

        assert_eq!(domain.len(), 3);
        assert_eq!(
            domain[1].as_f64(),
            Some(7.5),
            "threshold should replace domain[1]"
        );
    }
}

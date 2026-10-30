//! dp_analysis.rs
//!
//! Dynamic programming analysis for MSI estimation.
//! Computes data for pseudotime, distribution, and windowed heatmap outputs.
//!
//! This module provides:
//! 1. Data Structures (`RegionProbability`, `DpResult`, `AfEvolutionResult`,
//!    `WindowResult`, `FilteredRegions`, `AnalysisConfig`)
//! 2. Core DP Algorithm (`run_msi_dp`)
//! 3. Parallelization Setup (`setup_thread_pool`)
//! 4. Filtering Regions by AF (`filter_regions_by_af`)
//! 5. Distribution Primitives (`run_dp_for_regions`, `find_map_estimate`)
//! 6. MSI Metrics (`calculate_msi_metrics`)
//! 7. Window Utilities (`get_window_slice`)
//! 8. Windowed & Global Analysis (`run_windowed_analysis`, `run_global_analysis`)
//! 9. Orchestrator (`run_af_evolution_analysis`)

use std::collections::{HashMap, HashSet};

use anyhow::Result;
use rayon::prelude::*;

use super::extraction::{RegionSummary, Variant};
use super::OutputRequirements;
use crate::constants::EPSILON;
use crate::utils::stats::{calculate_percentage, usize_to_f64_exact};

/* ============ Data Structures =================== */

/// Region probability for DP computation.
///
/// Represents the probability that a microsatellite region is stable,
/// calculated as the product of all variant absence probabilities.
///
/// p_stable = P(all variants absent) = Π(prob_absent)
/// p_unstable = 1 - p_stable (computed in DP algorithm)
struct RegionProbability {
    /// P(all variants in the region absent): the product of their `prob_absent`.
    p_stable: f64,
}

/// DP result for one k value.
///
/// Represents the probability that exactly k regions are unstable,
/// along with the corresponding MSI score.
pub(super) struct DpResult {
    /// Number of unstable regions.
    pub(super) k: usize,
    /// MSI score for this k (k / total_regions × 100).
    pub(super) msi_score: f32,
    /// P(exactly k regions unstable).
    pub(super) probability: f64,
}

/// Result for one sample at one AF threshold.
///
/// Contains the MSI analysis results including MAP estimate,
/// uncertainty bounds, and optionally the full probability distribution.
pub(super) struct AfEvolutionResult {
    /// Total MS regions (denominator for the MSI score).
    pub(super) total_regions: usize,
    // Only computed if pseudotime output requested:
    /// MAP estimate of the number of unstable regions.
    pub(super) k_map: Option<usize>,
    /// MSI score at k_map (k_map / total_regions × 100).
    pub(super) msi_score_map: Option<f32>,
    /// Regions with at least one variant at or above the AF threshold.
    pub(super) regions_af_pass: Option<usize>,
    /// Lower bound (percent): (k_map - map_std_dev) / total_regions × 100, clamped at 0.
    pub(super) uncertainty_lower: Option<f32>,
    /// Upper bound (percent): (k_map + map_std_dev) / total_regions × 100, clamped at 100.
    pub(super) uncertainty_upper: Option<f32>,
    /// Standard deviation around the MAP, in region counts: sqrt(Σ (k - k_map)² × P(k)).
    /// Measured from k_map, not from the mean.
    pub(super) map_std_dev: Option<f64>,
    // Only computed if distribution output requested AND af_threshold == distribution_af:
    /// Full P(K=k) distribution.
    pub(super) distribution: Option<Vec<DpResult>>,
}

/// MSI score for one genomic window.
///
/// Forms one data point in the windowed plot:
/// - X-axis: genomic position (`window_start`)
/// - Y-axis: MSI score (`msi_score`)
/// - Color: posterior probability (`posterior_probability`)
///
/// Produced by `run_windowed_analysis` when a heatmap output is requested
/// (`--plot-heatmap` or `--data-heatmap`).
/// The fixed AF threshold used is configured via `--windowed-af`.
pub(super) struct WindowResult {
    /// Chromosome name.
    pub(super) chrom: String,
    /// Window start position (0-based, inclusive).
    pub(super) window_start: u64,
    /// Window end position (0-based, exclusive).
    pub(super) window_end: u64,
    /// MSI score for this window (k_map / regions_in_window × 100).
    pub(super) msi_score: f32,
    /// P(k=k_map) from DP distribution.
    pub(super) posterior_probability: f64,
    /// Total MS regions in this window - denominator for MSI score.
    pub(super) regions_in_window: usize,
}

/// Filtered regions view.
///
/// Uses a flat vector of variant references with region boundary indices
/// to avoid allocating thousands of small Vecs (~1-3 per region).
/// `filter_regions_by_af` leaves out regions with no variant passing the filter.
struct FilteredRegions<'a> {
    /// Flat vector of all variant references.
    ///
    /// ```text
    /// Example: [&V1, &V2, &V3, &V4, &V5]
    /// ```
    variants: Vec<&'a Variant>,

    /// Region boundaries (start indices).
    ///
    /// ```text
    /// Example: [0, 2, 3] means:
    ///   Region 0: variants[0..2] = [&V1, &V2]
    ///   Region 1: variants[2..3] = [&V3]
    ///   Region 2: variants[3..5] = [&V4, &V5]
    /// ```
    region_starts: Vec<usize>,
}

impl<'a> FilteredRegions<'a> {
    /// Get the variants for one region.
    ///
    /// # Arguments
    /// * `region_idx` - Index of the region, must be less than `len()`
    ///
    /// # Returns
    /// Slice of the variant references of region `region_idx`: from its start index
    /// up to, but not including, the next region's start index, or to the end of
    /// `variants` for the last region.
    ///
    /// # Panics
    /// Panics if `region_idx >= self.len()`, because `region_starts[region_idx]` is out of bounds.
    fn get_region(&self, region_idx: usize) -> &[&'a Variant] {
        let start = self.region_starts[region_idx];
        let end = self
            .region_starts
            .get(region_idx + 1)
            .copied()
            .unwrap_or(self.variants.len());
        &self.variants[start..end]
    }

    /// Number of regions in the view.
    ///
    /// # Returns
    /// Count of regions, which is the length of `region_starts`.
    fn len(&self) -> usize {
        self.region_starts.len()
    }
}

/// Holds everything needed to run one MSI calling analysis: how much data
/// is being analyzed, the sample name, the MSI-High classification
/// threshold, the AF thresholds to analyze, the two fixed AF values used
/// for the distribution and heatmap outputs, window size for heatmap analysis,
/// and the compute resources (thread count) to use.
pub(super) struct AnalysisConfig<'a> {
    /// Total number of MS regions (denominator for the pseudotime and distribution MSI scores).
    pub(super) total_regions: usize,
    /// Sample name being analyzed.
    pub(super) sample: &'a str,
    /// Threshold (percentage) at or above which a sample is classified MSI-High.
    pub(super) msi_high_threshold: f32,
    /// Allele frequency thresholds to compute pseudotime/distribution analysis data.
    pub(super) af_thresholds: Vec<f32>,
    /// Thread count for the rayon pool (`None` = rayon default).
    pub(super) num_threads: Option<usize>,
    /// Sliding window width (bp) for heatmap analysis.
    pub(super) window_size: u64,
    /// Fixed AF at which the full P(K=k) distribution is populated.
    pub(super) distribution_af: f32,
    /// Fixed AF used for windowed heatmap analysis.
    pub(super) windowed_af: f32,
}

/* =========== Core DP Algorithm ================== */

/// Execute dynamic programming algorithm for MSI variant probability distribution.
///
/// Implements matrix-based approach where:
/// - Matrix has n columns (regions 0 to n-1) and n+1 rows (unstable counts 0 to n)
/// - Our implementation keeps only the previous and current column, each a vector of length n+1
/// - The cell at column k and row i contains P(exactly i unstable regions among regions 0 to k)
///
/// # Algorithm
/// - **Initialize column 0**: `prev_col[0] = p_stable_0`, `prev_col[1] = p_unstable_0`
/// - **Recurrence**: `curr_col[i] = prev_col[i] × p_stable_k + prev_col[i-1] × p_unstable_k`
///
/// # Arguments
/// * `region_probs` - Regions with instability probabilities where:
///   - `p_stable` = P(all variants absent) = Π(prob_absent)
///   - `p_unstable` = 1 - `p_stable`
///
/// # Returns
/// - Probability distribution `[P(0 unstable), P(1 unstable), ..., P(n unstable)]`
///   where sum equals 1.0 (to rounding error) and each `P(i)` represents probability of exactly
///   i regions being unstable
/// - `[1.0]` if `region_probs` is empty
///
/// # Example
/// ```text
/// Input: 2 regions with p_unstable = [0.3, 0.4]
///
/// Step 1 (Region 0):
///   prev_col[0] = 0.7  // P(0 unstable) = p_stable_0
///   prev_col[1] = 0.3  // P(1 unstable) = p_unstable_0
///
/// Step 2 (Region 1):
///   p_stable_1 = 0.6, p_unstable_1 = 0.4
///   
///   curr_col[0] = prev_col[0] × 0.6 = 0.7 × 0.6 = 0.42
///   curr_col[1] = prev_col[1] × 0.6 + prev_col[0] × 0.4 = 0.3 × 0.6 + 0.7 × 0.4 = 0.46
///   curr_col[2] = prev_col[1] × 0.4 = 0.3 × 0.4 = 0.12
///
/// Result: [0.42, 0.46, 0.12]
/// Meaning: P(0 unstable)=42%, P(1 unstable)=46%, P(2 unstable)=12%
/// ```
///
/// # Note
/// - Uses space-optimized approach storing only previous column instead of full matrix
fn run_msi_dp(region_probs: &[RegionProbability]) -> Vec<f64> {
    let n = region_probs.len();

    // Base case: no regions means P(0 unstable) = 100%
    if n == 0 {
        return vec![1.0];
    }

    // Initialize probability vector: [P(0), P(1), P(2), ..., P(n)]
    let mut prev_col = vec![0.0; n + 1];

    // Column 0: Initialize with first region probabilities
    let p_stable_0 = region_probs[0].p_stable;
    let p_unstable_0 = 1.0 - p_stable_0;

    prev_col[0] = p_stable_0; // P(0 unstable) = first region stable
    prev_col[1] = p_unstable_0; // P(1 unstable) = first region unstable

    // Process remaining regions using recurrence relation
    #[allow(clippy::needless_range_loop)]
    for k in 1..n {
        let p_stable_k = region_probs[k].p_stable;
        let p_unstable_k = 1.0 - p_stable_k;

        let mut curr_col = vec![0.0; n + 1];

        // Base case: P(0 unstable) = previous P(0 unstable) × current region stable
        curr_col[0] = prev_col[0] * p_stable_k;

        // Recurrence: P(exactly i unstable) has two paths:
        // Path 1: Had i unstable, current region stable → still i unstable
        // Path 2: Had i-1 unstable, current region unstable → now i unstable
        for i in 1..=k + 1 {
            // Up to k+1 total unstable regions possible
            curr_col[i] = prev_col[i] * p_stable_k      // Path 1: Stay at i
                        + prev_col[i - 1] * p_unstable_k; // Path 2: Increment to i
        }

        prev_col = curr_col; // Update for next iteration
    }

    prev_col // Final distribution: [P(0), P(1), ..., P(n)]
}

/* =========== Parallelization Setup ============== */

/// Configure the rayon global thread pool.
///
/// Called once at the start of `run_af_evolution_analysis` before
/// any parallel work begins. Both `run_global_analysis` and
/// `run_windowed_analysis` share this pool automatically.
///
/// Keeps the existing pool, with a warning, if the global pool is already configured.
///
/// # Arguments
/// * `num_threads` - Thread count for the global pool, `None` leaves rayon's default
///
/// # Returns
/// * Nothing
fn setup_thread_pool(num_threads: Option<usize>) {
    if let Some(threads) = num_threads {
        // Note: In CLI usage, build_global() should always succeed on first call.
        // Error handling included for future library usage or testing scenarios.
        match rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
        {
            Ok(_) => {
                info!("Using {} threads (CLI specified)", threads);
            }
            Err(_) => {
                warn!(
                    "Using {} threads (global pool already configured, requested {} ignored)",
                    rayon::current_num_threads(),
                    threads
                );
            }
        }
    } else {
        info!(
            "Using {} threads (rayon default)",
            rayon::current_num_threads()
        );
    }
}

/* =========== Filtering Regions by AF ============ */

/// Filter a sample's regions down to variants meeting a minimum allele frequency.
///
/// Creates a flat filtered view avoiding Vec-per-region allocation overhead.
/// Optimized for MSI analysis where most regions have ~ 1-3 variants.
///
/// # Arguments
/// * `regions` - All regions with variants (lifetime 'a)
/// * `af_threshold` - Minimum AF
///
/// # Returns
/// Flat filtered view with region boundary tracking.
fn filter_regions_by_af<'a>(
    regions: &'a [RegionSummary],
    af_threshold: f32,
) -> FilteredRegions<'a> {
    let mut all_variants = Vec::new();
    let mut region_starts = Vec::new();

    for region in regions {
        let region_start = all_variants.len();
        let mut found_any = false;

        for variant in &region.variants {
            if variant.af >= af_threshold {
                all_variants.push(variant);
                found_any = true;
            }
        }

        if found_any {
            region_starts.push(region_start);
        }
    }

    FilteredRegions {
        variants: all_variants,
        region_starts,
    }
}

/* ============ Distribution Primitives =========== */

/// Compute DP probability distribution from pre-filtered regions.
///
/// Converts filtered regions into region stability probabilities
/// then runs the DP algorithm to get P(k regions unstable).
///
/// # Arguments
/// * `filtered` - Pre-filtered region view from `filter_regions_by_af`
///
/// # Returns
/// Probability distribution `[P(0 unstable), ..., P(n unstable)]`
/// summing to 1.0 (to rounding error). Returns `[1.0]` if no regions passed filtering.
fn run_dp_for_regions(filtered: &FilteredRegions) -> Vec<f64> {
    let region_probs: Vec<RegionProbability> = (0..filtered.len())
        .map(|i| {
            // P(all variants absent) = Π(p_absent)
            let p_all_absent: f64 = filtered
                .get_region(i)
                .iter()
                .map(|v| v.prob_absent)
                .product();

            RegionProbability {
                p_stable: p_all_absent,
            }
        })
        .collect();

    run_msi_dp(&region_probs)
}

/// Find MAP estimate, MSI score, and posterior probability from a DP distribution.
///
/// Finds k that maximizes P(K=k) - the maximum a posteriori estimate
/// of unstable region count. Ties broken by taking the highest k.
///
/// # Arguments
/// * `dist` - DP probability distribution from `run_dp_for_regions`
/// * `total_regions` - Denominator for the MSI percentage, either the genome-wide region
///   count (global analysis) or the region count within a single window (windowed analysis)
///
/// # Returns
/// `(k_map, msi_score, posterior_probability)` where:
/// - `k_map` - the k maximizing P(K=k), the highest such k on ties
/// - `msi_score = k_map / total_regions × 100`
/// - `posterior_probability = P(k=k_map)` - confidence in the estimate
///
/// # Panics
/// Panics if `dist` is empty, if `dist` has two or more entries and one is NaN, or if
/// `total_regions` is 0.
fn find_map_estimate(dist: &[f64], total_regions: usize) -> (usize, f32, f64) {
    let k_map = dist
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(k, _)| k)
        .expect("dist is guaranteed non-empty: run_dp_for_regions always returns at least [1.0]");

    let msi_score = calculate_percentage(k_map, total_regions);
    let posterior_probability = dist[k_map];

    (k_map, msi_score, posterior_probability)
}

/* ============ MSI Metrics ======================= */

/// Calculate MSI metrics for filtered regions.
///
/// Computes region instability probabilities, runs DP, calculates MAP estimate,
/// and computes uncertainty bounds if requested.
///
/// # Steps
/// 1. **Region Probabilities**: P(region unstable) = 1 - Π(p_absent)
/// 2. **DP**: Compute probability distribution P(K = k)
/// 3. **MAP** (if pseudotime requested): Find k that maximizes P(K = k)
/// 4. **MSI Score** (if pseudotime requested): (k_map / total_regions) × 100
/// 5. **Uncertainty** (if pseudotime requested): Standard deviation bounds
/// 6. **Distribution** (if requested and `af_threshold == distribution_af`): Full P(K = k) array
///
/// # Arguments
/// * `filtered` - Filtered view of variants passing AF threshold
/// * `total_regions` - Total number of MS regions
/// * `af_threshold` - AF threshold used for filtering
/// * `output_req` - Which outputs are requested (controls what to compute)
/// * `distribution_af` - Fixed AF at which `distribution` is populated,
///   `distribution` is only set when `output_req.needs_distribution` is true
///   and `af_threshold` equals this value
///
/// # Returns
/// `AfEvolutionResult` for this AF threshold, where:
/// - `total_regions` - always set
/// - `k_map`, `msi_score_map`, `regions_af_pass`, `uncertainty_lower`, `uncertainty_upper`
///   and `map_std_dev` - set only if `output_req.needs_pseudotime`, otherwise `None`
/// - `distribution` - set only if `output_req.needs_distribution` and `af_threshold`
///   equals `distribution_af`, otherwise `None`
///
/// # Panics
/// Panics if `total_regions` is 0 and a percentage is computed (pseudotime requested,
/// or the distribution is populated). Unreachable from `call_msi`, which returns
/// `MsiBedRegionsEmpty` before this point.
fn calculate_msi_metrics(
    filtered: &FilteredRegions,
    total_regions: usize,
    af_threshold: f32,
    output_req: OutputRequirements,
    distribution_af: f32,
) -> AfEvolutionResult {
    // Steps 1-2: Per-region stability probability and the DP distribution,
    // both computed through run_dp_for_regions
    let distribution_raw = run_dp_for_regions(filtered);

    // Steps 3-5: MAP estimate, MSI score and uncertainty bounds, only if pseudotime is requested
    let (k_map, msi_score_map, regions_af_pass, uncertainty_lower, uncertainty_upper, map_std_dev) =
        if output_req.needs_pseudotime {
            let (k_map_raw, msi_score_map_raw, _) =
                find_map_estimate(&distribution_raw, total_regions);
            let regions_af_pass_count = filtered.len();

            let k_map_f64 = usize_to_f64_exact(k_map_raw);

            // Mean squared deviation from k_map (not from the mean): Σ[(k - k_map)² × P(k)]
            let variance_f64: f64 = distribution_raw
                .iter()
                .enumerate()
                .map(|(k, &prob)| {
                    let k_f64 = usize_to_f64_exact(k);
                    let diff = k_f64 - k_map_f64;
                    let diff_squared = diff * diff;
                    diff_squared * prob
                })
                .sum();

            // Standard deviation around k_map
            let std_dev_f64 = variance_f64.sqrt();

            // Bounds k_map ± std dev, clamped to [0, total_regions], as percentages
            let total_f64 = usize_to_f64_exact(total_regions);
            let lower_k = (k_map_f64 - std_dev_f64).max(0.0);
            let upper_k = (k_map_f64 + std_dev_f64).min(total_f64);
            let lower = ((lower_k / total_f64) * 100.0) as f32;
            let upper = ((upper_k / total_f64) * 100.0) as f32;

            (
                Some(k_map_raw),
                Some(msi_score_map_raw),
                Some(regions_af_pass_count),
                Some(lower),
                Some(upper),
                Some(std_dev_f64),
            )
        } else {
            (None, None, None, None, None, None)
        };

    // Step 6: Create full distribution (only if distribution output requested
    // and current AF threshold matches distribution_af)
    let distribution = if output_req.needs_distribution
        && (af_threshold - distribution_af).abs() < EPSILON as f32
    {
        Some(
            distribution_raw
                .iter()
                .enumerate()
                .map(|(k, &probability)| {
                    let msi_score = calculate_percentage(k, total_regions);
                    DpResult {
                        k,
                        msi_score,
                        probability,
                    }
                })
                .collect(),
        )
    } else {
        None
    };

    AfEvolutionResult {
        total_regions,
        k_map,
        msi_score_map,
        regions_af_pass,
        uncertainty_lower,
        uncertainty_upper,
        map_std_dev,
        distribution,
    }
}

/* ============ Window Utilities ================== */

/// Get a zero-copy sub-slice of regions falling within a genomic window.
///
/// # Precondition
/// `regions` must be sorted by `(chrom, start)`. This should be guaranteed by
/// `extract_regions`, otherwise results would be incorrect.
///
/// # Arguments
/// * `regions` - All extracted regions sorted by `(chrom, start)`
/// * `chrom` - Chromosome to filter for
/// * `window_start` - Window start position (0-based, inclusive)
/// * `window_end` - Window end position (0-based, exclusive)
///
/// # Returns
/// Sub-slice of `regions` whose chrom matches and start falls in
/// `[window_start, window_end)`. Returns empty slice if no match.
fn get_window_slice<'a>(
    regions: &'a [RegionSummary],
    chrom: &str,
    window_start: u64,
    window_end: u64,
) -> &'a [RegionSummary] {
    let start_idx =
        regions.partition_point(|r| (r.chrom.as_str(), r.start) < (chrom, window_start));
    let end_idx = start_idx
        + regions[start_idx..]
            .partition_point(|r| (r.chrom.as_str(), r.start) < (chrom, window_end));

    &regions[start_idx..end_idx]
}

/* ============ Windowed & Global Analysis ======== */

/// Compute windowed MSI scores across the genome.
///
/// Splits each chromosome into consecutive, non-overlapping windows of
/// `window_size` bases and computes MSI score and posterior probability per
/// window using a single fixed AF threshold. A region belongs to the window
/// containing its start position. Only windows containing at least one MS
/// region are included.
///
/// Per-window denominator (`regions_in_window`) makes windows
/// spatially comparable - each score reflects instability fraction
/// within that window, not genome-wide.
///
/// # Arguments
/// * `regions` - All extracted regions sorted by `(chrom, start)`
/// * `af_threshold` - Fixed AF threshold from `--windowed-af`
/// * `window_size` - Window width in bases (e.g. 1_000_000 for 1Mb)
///
/// # Returns
/// Vec of `WindowResult` ordered by `(chrom, window_start)`, with chromosome
/// names compared as strings (`chr1`, `chr10`, `chr2`).
///
/// # Panics
/// In debug builds, panics (failed `debug_assert!`) if `window_size` is 0 or
/// `regions` is empty.
fn run_windowed_analysis(
    regions: &[RegionSummary],
    af_threshold: f32,
    window_size: u64,
) -> Vec<WindowResult> {
    debug_assert!(
        window_size != 0,
        "window_size=0 guaranteed unreachable: MSIConfig::validate() requires \
        sliding_window != 0 whenever a heatmap output is requested, and this \
        function is only called when needs_heatmap is true"
    );

    debug_assert!(
        !regions.is_empty(),
        "regions=[] guaranteed unreachable: call_msi returns MsiBedRegionsEmpty before \
        reaching run_windowed_analysis whenever total_ms_regions == 0, and this function \
        is only called when needs_heatmap is true"
    );

    // Collect unique chromosomes
    let mut chroms: HashSet<&str> = HashSet::new();
    for region in regions {
        chroms.insert(region.chrom.as_str());
    }

    // Build (chrom, window_start) work items - single AF
    let work_items: Vec<(&str, u64)> = chroms
        .iter()
        .flat_map(|&chrom| {
            let max_start = regions
                .iter()
                .filter(|r| r.chrom == chrom)
                .map(|r| r.start)
                .max()
                .unwrap_or(0);
            let highest_window_idx = max_start / window_size;
            let n_windows = highest_window_idx + 1;
            (0..n_windows).map(move |w| (chrom, w * window_size))
        })
        .collect();

    info!(
        "Running windowed analysis: {} windows across {} chromosomes",
        work_items.len(),
        chroms.len()
    );

    let mut all_results: Vec<WindowResult> = work_items
        .par_iter()
        .filter_map(|&(chrom, window_start)| {
            let window_end = window_start + window_size;

            // Zero-copy sub-slice for this window
            let window_slice = get_window_slice(regions, chrom, window_start, window_end);
            if window_slice.is_empty() {
                return None;
            }

            let regions_in_window = window_slice.len();

            let filtered = filter_regions_by_af(window_slice, af_threshold);

            let dist = run_dp_for_regions(&filtered);
            let (_, msi_score, posterior_probability) = find_map_estimate(&dist, regions_in_window);

            Some(WindowResult {
                chrom: chrom.to_string(),
                window_start,
                window_end,
                msi_score,
                posterior_probability,
                regions_in_window,
            })
        })
        .collect();

    // Sort by (chrom, window_start)
    all_results.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.window_start.cmp(&b.window_start))
    });

    info!(
        "Windowed analysis complete: {} non-empty windows",
        all_results.len()
    );

    all_results
}

/// Run parallel MSI analysis across all AF thresholds.
/// Responsible for generating data for pseudotime and distribution outputs.
///
/// # Arguments
/// * `regions`            - Regions with variants from extraction
/// * `total_regions`      - Total MS regions (denominator for MSI score)
/// * `af_thresholds`      - AF thresholds to analyze in parallel
/// * `distribution_af`    - Fixed AF forwarded to `calculate_msi_metrics` on each call
///   marking which AF should populate the full distribution conditionally
/// * `output_req`         - Which outputs to compute
///
/// # Returns
/// HashMap keyed by AF threshold string (e.g. "0.00"), MSI result.
fn run_global_analysis(
    regions: &[RegionSummary],
    total_regions: usize,
    af_thresholds: &[f32],
    distribution_af: f32,
    output_req: OutputRequirements,
) -> HashMap<String, AfEvolutionResult> {
    info!(
        "Running global analysis: {} AF thresholds in parallel",
        af_thresholds.len()
    );

    let final_results: HashMap<String, AfEvolutionResult> = af_thresholds
        .par_iter()
        .map(|af_threshold| {
            let filtered = filter_regions_by_af(regions, *af_threshold);
            let result = calculate_msi_metrics(
                &filtered,
                total_regions,
                *af_threshold,
                output_req,
                distribution_af,
            );

            (format!("{:.2}", af_threshold), result)
        })
        .collect();

    info!("Global analysis complete: {} results", final_results.len());
    final_results
}

/* ============ Orchestrator ====================== */

/// Run AF evolution analysis across AF thresholds for one sample.
///
/// Performs parallel MSI analysis across the AF thresholds (for
/// pseudotime/distribution outputs) and, if requested, windowed heatmap
/// analysis across the genome.
///
/// # Arguments
/// * `regions`    - Regions with variants from extraction
/// * `config`     - Config items: total regions of interest, sample, msi threshold,
///   AF Thresholds for pseudotime analysis, fixed distribution/windowed AF values,
///   windowed heatmap window size and compute resource
/// * `output_req` - Which outputs requested (controls computation)
///
/// # Note
/// This function does not trim `config.af_thresholds` itself - it computes
/// the DP for exactly whatever thresholds are present. The caller
/// (`call_msi` in `mod.rs`) is responsible for pre-trimming the list to
/// avoid wasted computation: the full list when pseudotime is requested,
/// a single-element `[distribution_af]` when only distribution is
/// requested, or empty when neither is needed.
///
/// # Returns
/// Tuple of `(global_results, window_results)` where:
/// - `global_results`: HashMap keyed by AF threshold string -> MSI result.
///   Empty if neither pseudotime nor distribution requested.
/// - `window_results`: `Vec<WindowResult>`, empty if heatmap not requested
///   or no windows contained regions.
///
/// # Example
/// assert_eq!(run_af_evolution_analysis(&regions, config, output_req), (global_results, window_results));
pub(super) fn run_af_evolution_analysis(
    regions: &[RegionSummary],
    config: AnalysisConfig<'_>,
    output_req: OutputRequirements,
) -> Result<(HashMap<String, AfEvolutionResult>, Vec<WindowResult>)> {
    info!("Sample: {:?}", config.sample);
    info!("AF thresholds: {:?}", config.af_thresholds);
    info!("MSI-High threshold: {}%", config.msi_high_threshold);
    info!("Total regions (BED): {}", config.total_regions);
    info!("Regions with variants: {}", regions.len());
    info!("Output requirements:");
    info!("    - Pseudotime: {}", output_req.needs_pseudotime);
    if output_req.needs_distribution {
        info!(
            "    - Distribution (AF={}): {}",
            config.distribution_af, output_req.needs_distribution
        );
    } else {
        info!("    - Distribution: {}", output_req.needs_distribution);
    }
    if output_req.needs_heatmap {
        info!(
            "    - Heatmap (AF={}): {}",
            config.windowed_af, output_req.needs_heatmap
        );
    } else {
        info!("    - Heatmap: {}", output_req.needs_heatmap);
    }

    // Step 1: Configure rayon thread pool
    setup_thread_pool(config.num_threads); // once shared by both analysis functions

    // Step 2: Run global analysis — skipped if only heatmap requested
    let global_results = if output_req.needs_pseudotime || output_req.needs_distribution {
        run_global_analysis(
            regions,
            config.total_regions,
            &config.af_thresholds,
            config.distribution_af,
            output_req,
        )
    } else {
        HashMap::new()
    };

    // Step 3: Run windowed analysis — only if heatmap requested
    let window_results = if output_req.needs_heatmap {
        run_windowed_analysis(regions, config.windowed_af, config.window_size)
    } else {
        Vec::new()
    };

    Ok((global_results, window_results))
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::constants::test_constants::{TEST_EPSILON, TEST_EPSILON_F32};

    /// Helper: Create a test Variant with given absence probability and allele frequency.
    fn make_variant(prob: f64, af: f32) -> Variant {
        Variant {
            prob_absent: prob,
            af,
        }
    }

    /// Helper: Create a test RegionSummary with full control over all fields.
    fn make_region(
        chrom: &str,
        start: u64,
        variants: Vec<Variant>,
        has_real_indel: bool,
    ) -> RegionSummary {
        RegionSummary {
            chrom: chrom.to_string(),
            start,
            variants,
            has_real_indel,
        }
    }

    /// Helper: Create a test RegionSummary with default location fields,
    /// for tests that only care about variants, not genomic position.
    fn make_region_simple(variants: Vec<Variant>) -> RegionSummary {
        make_region("chr1", 0, variants, true)
    }

    /* ============ FilteredRegions tests ============ */

    #[test]
    fn test_filtered_regions_get_region() {
        let v1 = make_variant(0.01, 0.8);
        let v2 = make_variant(0.02, 0.9);
        let v3 = make_variant(0.03, 0.7);

        let filtered = FilteredRegions {
            variants: vec![&v1, &v2, &v3],
            region_starts: vec![0, 2],
        };

        // Region 0: variants[0..2]
        let r0 = filtered.get_region(0);
        assert_eq!(r0.len(), 2);
        assert!((r0[0].prob_absent - 0.01).abs() < TEST_EPSILON);
        assert!((r0[1].prob_absent - 0.02).abs() < TEST_EPSILON);

        // Region 1: variants[2..3]
        let r1 = filtered.get_region(1);
        assert_eq!(r1.len(), 1);
        assert!((r1[0].prob_absent - 0.03).abs() < TEST_EPSILON);
    }

    #[test]
    #[should_panic(expected = "index out of bounds: the len is 2 but the index is 2")]
    fn test_filtered_regions_get_region_out_of_range() {
        let v1 = make_variant(0.01, 0.8);
        let v2 = make_variant(0.02, 0.9);

        let filtered = FilteredRegions {
            variants: vec![&v1, &v2],
            region_starts: vec![0, 1],
        };

        // 2 regions, so valid indices are 0 and 1: index 2 == len() is out of range
        filtered.get_region(2);
    }

    #[test]
    fn test_filtered_regions_len() {
        let v1 = make_variant(0.01, 0.8);

        let filtered = FilteredRegions {
            variants: vec![&v1],
            region_starts: vec![0],
        };

        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered.get_region(0).len(), 1);
    }

    #[test]
    fn test_filtered_regions_len_empty() {
        let filtered = FilteredRegions {
            variants: vec![],
            region_starts: vec![],
        };
        assert_eq!(filtered.len(), 0);
    }

    /* ============== run_msi_dp tests =============== */

    #[test]
    fn test_run_msi_dp_empty() {
        let probs = vec![];
        let dist = run_msi_dp(&probs);

        assert_eq!(dist.len(), 1);
        assert!((dist[0] - 1.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_msi_dp_single_region() {
        let probs = vec![RegionProbability { p_stable: 0.7 }];
        let dist = run_msi_dp(&probs);

        assert_eq!(dist.len(), 2);
        assert!((dist[0] - 0.7).abs() < TEST_EPSILON); // P(0 unstable) = 0.7
        assert!((dist[1] - 0.3).abs() < TEST_EPSILON); // P(1 unstable) = 0.3
    }

    #[test]
    fn test_run_msi_dp_two_regions() {
        let probs = vec![
            RegionProbability { p_stable: 0.7 },
            RegionProbability { p_stable: 0.6 },
        ];
        let dist = run_msi_dp(&probs);

        assert_eq!(dist.len(), 3);
        // P(0) = 0.7 × 0.6 = 0.42
        // P(1) = 0.7 × 0.4 + 0.3 × 0.6 = 0.46
        // P(2) = 0.3 × 0.4 = 0.12
        assert!((dist[0] - 0.42).abs() < TEST_EPSILON);
        assert!((dist[1] - 0.46).abs() < TEST_EPSILON);
        assert!((dist[2] - 0.12).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_msi_dp_three_regions() {
        let probs = vec![
            RegionProbability { p_stable: 0.1 },
            RegionProbability { p_stable: 0.5 },
            RegionProbability { p_stable: 0.8 },
        ];
        let dist = run_msi_dp(&probs);

        assert_eq!(dist.len(), 4);
        assert!((dist[0] - 0.04).abs() < TEST_EPSILON);
        assert!((dist[1] - 0.41).abs() < TEST_EPSILON);
        assert!((dist[2] - 0.46).abs() < TEST_EPSILON);
        assert!((dist[3] - 0.09).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_msi_dp_distribution_sums_to_one() {
        let probs = vec![
            RegionProbability { p_stable: 0.2 },
            RegionProbability { p_stable: 0.5 },
            RegionProbability { p_stable: 0.8 },
        ];
        let dist = run_msi_dp(&probs);

        let sum: f64 = dist.iter().sum();
        assert!((sum - 1.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_msi_dp_all_stable() {
        let probs = vec![
            RegionProbability { p_stable: 1.0 },
            RegionProbability { p_stable: 1.0 },
        ];
        let dist = run_msi_dp(&probs);
        assert!((dist[0] - 1.0).abs() < TEST_EPSILON);
        assert!((dist[1] - 0.0).abs() < TEST_EPSILON);
        assert!((dist[2] - 0.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_msi_dp_all_unstable() {
        let probs = vec![
            RegionProbability { p_stable: 0.0 },
            RegionProbability { p_stable: 0.0 },
        ];
        let dist = run_msi_dp(&probs);
        assert!((dist[0] - 0.0).abs() < TEST_EPSILON);
        assert!((dist[1] - 0.0).abs() < TEST_EPSILON);
        assert!((dist[2] - 1.0).abs() < TEST_EPSILON);
    }

    /* ========== setup_thread_pool tests ============ */

    #[test]
    fn test_setup_thread_pool_none_does_not_panic() {
        setup_thread_pool(None);
    }

    #[test]
    fn test_setup_thread_pool_already_configured_keeps_pool() {
        // Run a task on the global pool so it certainly exists:
        // build_global() fails once the pool exists
        rayon::spawn(|| {});
        let before = rayon::current_num_threads();
        setup_thread_pool(Some(before + 1));
        assert_eq!(rayon::current_num_threads(), before);
    }

    #[test]
    fn test_setup_thread_pool_pool_created_by_thread_count_read_keeps_pool() {
        let before = rayon::current_num_threads();
        setup_thread_pool(Some(before + 1));
        assert_eq!(rayon::current_num_threads(), before);
    }

    /* ========= filter_regions_by_af tests ========== */

    #[test]
    fn test_filter_regions_by_af_basic() {
        let regions = vec![make_region_simple(vec![
            make_variant(0.01, 0.8), // Pass
            make_variant(0.02, 0.4), // Fail
            make_variant(0.03, 0.7), // Pass
        ])];

        let filtered = filter_regions_by_af(&regions, 0.6);

        assert_eq!(filtered.len(), 1);
        let region_1 = filtered.get_region(0);
        assert_eq!(region_1.len(), 2); // V1 and V3
        assert!((region_1[0].prob_absent - 0.01).abs() < TEST_EPSILON);
        assert!((region_1[1].prob_absent - 0.03).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_filter_regions_by_af_zero_includes_all() {
        let regions = vec![make_region_simple(vec![
            make_variant(0.01, 0.0),
            make_variant(0.02, 0.5),
            make_variant(0.03, 1.0),
        ])];

        let filtered = filter_regions_by_af(&regions, 0.0);

        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered.get_region(0).len(), 3);
    }

    #[test]
    fn test_filter_regions_by_af_multiple_regions_some_filtered_out() {
        let regions = vec![
            make_region("chr1", 100, vec![make_variant(0.01, 0.8)], true),
            make_region("chr1", 200, vec![make_variant(0.02, 0.3)], true), // Filtered out
            make_region("chr1", 300, vec![make_variant(0.03, 0.9)], true),
        ];

        let filtered = filter_regions_by_af(&regions, 0.5);
        assert_eq!(filtered.len(), 2); // Only first and third
        assert!((filtered.get_region(0)[0].af - 0.8).abs() < TEST_EPSILON_F32);
        assert!((filtered.get_region(1)[0].af - 0.9).abs() < TEST_EPSILON_F32);
    }

    #[test]
    fn test_filter_regions_by_af_boundary_exact_threshold() {
        // af == threshold is kept: >= is inclusive, not strict >
        let regions = vec![make_region_simple(vec![make_variant(0.01, 0.6)])];
        let filtered = filter_regions_by_af(&regions, 0.6);
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered.get_region(0).len(), 1);
    }

    #[test]
    fn test_filter_regions_by_af_empty_regions_input() {
        let regions: Vec<RegionSummary> = vec![];
        let filtered = filter_regions_by_af(&regions, 0.5);
        assert_eq!(filtered.len(), 0);
    }

    /* ========== run_dp_for_regions tests =========== */

    #[test]
    fn test_run_dp_for_regions_empty() {
        let filtered = FilteredRegions {
            variants: vec![],
            region_starts: vec![],
        };
        let dist = run_dp_for_regions(&filtered);
        assert_eq!(dist, vec![1.0]);
    }

    #[test]
    fn test_run_dp_for_regions_single_region() {
        // prob_absent=0.7, p_stable=0.7
        let v = make_variant(0.7, 0.8);
        let filtered = FilteredRegions {
            variants: vec![&v],
            region_starts: vec![0],
        };
        let dist = run_dp_for_regions(&filtered);
        assert_eq!(dist.len(), 2);
        assert!((dist[0] - 0.7).abs() < TEST_EPSILON); // P(0 unstable)
        assert!((dist[1] - 0.3).abs() < TEST_EPSILON); // P(1 unstable)
    }

    #[test]
    fn test_run_dp_for_regions_sums_to_one() {
        let v1 = make_variant(0.3, 0.8);
        let v2 = make_variant(0.5, 0.9);
        let v3 = make_variant(0.7, 0.7);
        let filtered = FilteredRegions {
            variants: vec![&v1, &v2, &v3],
            region_starts: vec![0, 1, 2],
        };
        let dist = run_dp_for_regions(&filtered);
        let sum: f64 = dist.iter().sum();
        assert!((sum - 1.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_dp_for_regions_multiple_variants_per_region() {
        // Region with two variants: p_stable = 0.3 × 0.5 = 0.15
        let v1 = make_variant(0.3, 0.8);
        let v2 = make_variant(0.5, 0.9);
        let filtered = FilteredRegions {
            variants: vec![&v1, &v2],
            region_starts: vec![0], // one region with two variants
        };
        let dist = run_dp_for_regions(&filtered);
        assert_eq!(dist.len(), 2);
        assert!((dist[0] - 0.15).abs() < TEST_EPSILON); // P(0) = 0.15
        assert!((dist[1] - 0.85).abs() < TEST_EPSILON); // P(1) = 0.85
    }

    #[test]
    fn test_run_dp_for_regions_uneven_region_sizes() {
        // Region 0 = [v1, v2]: p_stable = 0.5 × 0.4 = 0.2
        // Region 1 = [v3]:     p_stable = 0.7
        let v1 = make_variant(0.5, 0.8);
        let v2 = make_variant(0.4, 0.9);
        let v3 = make_variant(0.7, 0.7);
        let filtered = FilteredRegions {
            variants: vec![&v1, &v2, &v3],
            region_starts: vec![0, 2],
        };
        let dist = run_dp_for_regions(&filtered);

        assert_eq!(dist.len(), 3);
        // P(0) = 0.2 × 0.7 = 0.14
        // P(1) = 0.8 × 0.7 + 0.2 × 0.3 = 0.56 + 0.06 = 0.62
        // P(2) = 0.8 × 0.3 = 0.24
        assert!((dist[0] - 0.14).abs() < TEST_EPSILON);
        assert!((dist[1] - 0.62).abs() < TEST_EPSILON);
        assert!((dist[2] - 0.24).abs() < TEST_EPSILON);
    }

    /* ========== find_map_estimate tests ============ */

    #[test]
    fn test_find_map_estimate_basic() {
        // P(0)=0.42, P(1)=0.46, P(2)=0.12, k_map=1
        let dist = vec![0.42, 0.46, 0.12];
        let (k_map, score, posterior_probability) = find_map_estimate(&dist, 2);
        assert_eq!(k_map, 1);
        assert!((score - 50.0).abs() < TEST_EPSILON_F32); // 1/2 × 100
        assert!((posterior_probability - 0.46).abs() < TEST_EPSILON); // P(k=1)
    }

    #[test]
    fn test_find_map_estimate_tie_takes_higher_k() {
        // Equal probabilities, higher k wins
        let dist = vec![0.5, 0.5];
        let (k_map, _, posterior_probability) = find_map_estimate(&dist, 10);
        assert_eq!(k_map, 1); // higher k wins on tie
        assert!((posterior_probability - 0.5).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_find_map_estimate_all_stable() {
        // P(0)=1.0, k_map=0, score=0%
        let dist = vec![1.0, 0.0, 0.0];
        let (k_map, score, posterior_probability) = find_map_estimate(&dist, 2);
        assert_eq!(k_map, 0);
        assert_eq!(score, 0.0);
        assert!((posterior_probability - 1.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_find_map_estimate_all_unstable() {
        // P(2)=1.0, k_map=2, score=100%
        let dist = vec![0.0, 0.0, 1.0];
        let (k_map, score, posterior_probability) = find_map_estimate(&dist, 2);
        assert_eq!(k_map, 2);
        assert!((score - 100.0).abs() < TEST_EPSILON_F32);
        assert!((posterior_probability - 1.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_find_map_estimate_total_regions_differs_from_dist_len() {
        // k_map=2, total=10, 2/10 × 100 = 20%
        // posterior_probability is dist[k_map] = 0.7, the largest entry, so no separate test is needed
        let dist = vec![0.1, 0.2, 0.7];
        let (k_map, score, posterior_probability) = find_map_estimate(&dist, 10);
        assert_eq!(k_map, 2);
        assert!((score - 20.0).abs() < TEST_EPSILON_F32);
        assert!((posterior_probability - 0.7).abs() < TEST_EPSILON);
    }

    #[test]
    #[should_panic(expected = "dist is guaranteed non-empty")]
    fn test_find_map_estimate_empty_dist() {
        let dist: Vec<f64> = vec![];
        find_map_estimate(&dist, 10);
    }

    #[test]
    #[should_panic(expected = "called `Option::unwrap()` on a `None` value")]
    fn test_find_map_estimate_nan_in_dist() {
        // partial_cmp returns None when NaN is compared with another value
        let dist = vec![f64::NAN, 0.5];
        find_map_estimate(&dist, 10);
    }

    #[test]
    #[should_panic(expected = "denominator must not be zero")]
    fn test_find_map_estimate_zero_total_regions() {
        let dist = vec![0.5, 0.5];
        find_map_estimate(&dist, 0);
    }

    /* ======== calculate_msi_metrics tests ========== */

    #[test]
    fn test_calculate_msi_metrics_no_outputs_requested() {
        let filtered = FilteredRegions {
            variants: vec![],
            region_starts: vec![],
        };

        // Unreachable from the CLI: MSIConfig::validate() returns
        // Error::MsiConfigOutputMissing when no output is requested.
        // Kept as a guard: nothing requested leaves the fields None.
        let output_req = OutputRequirements {
            needs_pseudotime: false,
            needs_distribution: false,
            needs_heatmap: false,
        };

        // af_threshold equals distribution_af, but nothing is requested
        let result = calculate_msi_metrics(&filtered, 100, 0.0, output_req, 0.0);

        assert!(result.k_map.is_none());
        assert!(result.msi_score_map.is_none());
        assert!(result.regions_af_pass.is_none());
        assert!(result.uncertainty_lower.is_none());
        assert!(result.uncertainty_upper.is_none());
        assert!(result.map_std_dev.is_none());
        assert!(result.distribution.is_none());
    }

    #[test]
    fn test_calculate_msi_metrics_empty_filtered_set() {
        let filtered = FilteredRegions {
            variants: vec![],
            region_starts: vec![],
        };
        let output_req = OutputRequirements {
            needs_pseudotime: true,
            needs_distribution: true,
            needs_heatmap: false,
        };

        let result = calculate_msi_metrics(&filtered, 100, 0.5, output_req, 0.5);

        // Empty filtered view: the DP distribution is [1.0],
        // so k_map = 0 and the distribution has one entry
        assert_eq!(result.k_map, Some(0));
        assert_eq!(result.regions_af_pass, Some(0));
        assert_eq!(result.msi_score_map, Some(0.0));
        assert_eq!(result.map_std_dev, Some(0.0));
        assert_eq!(result.uncertainty_lower, Some(0.0));
        assert_eq!(result.uncertainty_upper, Some(0.0));

        let dist = result.distribution.unwrap();
        assert_eq!(dist.len(), 1);
        assert_eq!(dist[0].k, 0);
        assert_eq!(dist[0].msi_score, 0.0);
        assert_eq!(dist[0].probability, 1.0);
    }

    #[test]
    fn test_calculate_msi_metrics_with_uncertainty() {
        let v1 = make_variant(0.1, 0.8);
        let v2 = make_variant(0.2, 0.9);

        let filtered = FilteredRegions {
            variants: vec![&v1, &v2],
            region_starts: vec![0, 1],
        };

        let output_req = OutputRequirements {
            needs_pseudotime: true,
            needs_distribution: false,
            needs_heatmap: false,
        };

        let result = calculate_msi_metrics(&filtered, 100, 0.5, output_req, 0.05);

        // p_stable = [0.1, 0.2]: P(0) = 0.1 × 0.2 = 0.02,
        // P(1) = 0.9 × 0.2 + 0.1 × 0.8 = 0.26, P(2) = 0.9 × 0.8 = 0.72
        // k_map = 2 and std around k_map = sqrt(2² × 0.02 + 1² × 0.26) = sqrt(0.34)
        let std_expected = (4.0 * 0.02 + 0.26_f64).sqrt();

        assert_eq!(result.k_map, Some(2));
        assert_eq!(result.regions_af_pass, Some(2));
        assert!((result.msi_score_map.unwrap() - 2.0).abs() < TEST_EPSILON_F32); // 2/100 × 100
        assert!((result.map_std_dev.unwrap() - std_expected).abs() < TEST_EPSILON);
        assert!(
            (result.uncertainty_lower.unwrap() - (2.0 - std_expected) as f32).abs()
                < TEST_EPSILON_F32
        );
        assert!(
            (result.uncertainty_upper.unwrap() - (2.0 + std_expected) as f32).abs()
                < TEST_EPSILON_F32
        );
        assert!(result.distribution.is_none());
    }

    #[test]
    fn test_calculate_msi_metrics_distribution_only_at_af_zero() {
        let v1 = make_variant(0.1, 0.8);

        let filtered = FilteredRegions {
            variants: vec![&v1],
            region_starts: vec![0],
        };

        let output_req = OutputRequirements {
            needs_pseudotime: false,
            needs_distribution: true,
            needs_heatmap: false,
        };

        // AF=0.0 should include distribution
        let result_af0 = calculate_msi_metrics(&filtered, 100, 0.0, output_req, 0.00);
        assert!(result_af0.distribution.is_some());
        let dist = result_af0.distribution.unwrap();
        assert_eq!(dist.len(), 2);
        let prob_sum: f64 = dist.iter().map(|d| d.probability).sum();
        assert!((prob_sum - 1.0).abs() < TEST_EPSILON);
        assert_eq!(dist[0].k, 0);
        assert_eq!(dist[1].k, 1);
        assert_eq!(dist[0].msi_score, 0.0); // 0/100 × 100
        assert!((dist[1].msi_score - 1.0).abs() < TEST_EPSILON_F32); // 1/100 × 100
        assert!((dist[0].probability - 0.1).abs() < TEST_EPSILON);
        assert!((dist[1].probability - 0.9).abs() < TEST_EPSILON);

        // AF=0.5 should NOT include distribution
        let result_af5 = calculate_msi_metrics(&filtered, 100, 0.5, output_req, 0.00);
        assert!(result_af5.distribution.is_none());
    }

    #[test]
    #[should_panic(expected = "denominator must not be zero")]
    fn test_calculate_msi_metrics_zero_total_regions() {
        let filtered = FilteredRegions {
            variants: vec![],
            region_starts: vec![],
        };
        let output_req = OutputRequirements {
            needs_pseudotime: true,
            needs_distribution: false,
            needs_heatmap: false,
        };

        calculate_msi_metrics(&filtered, 0, 0.0, output_req, 0.05);
    }

    #[test]
    fn test_calculate_msi_metrics_lower_bound_clamped_at_zero() {
        let v1 = make_variant(0.9, 0.8);
        let v2 = make_variant(0.9, 0.8);
        let filtered = FilteredRegions {
            variants: vec![&v1, &v2],
            region_starts: vec![0, 1],
        };
        let output_req = OutputRequirements {
            needs_pseudotime: true,
            needs_distribution: false,
            needs_heatmap: false,
        };

        let result = calculate_msi_metrics(&filtered, 2, 0.5, output_req, 0.05);

        // p_stable = [0.9, 0.9]: P(0) = 0.81, P(1) = 0.18, P(2) = 0.01, so k_map = 0
        // std around k_map = sqrt(1² × 0.18 + 2² × 0.01) = sqrt(0.22)
        let std_expected = (0.18 + 4.0 * 0.01_f64).sqrt();

        assert_eq!(result.k_map, Some(0));
        // k_map - std < 0, so the lower bound is clamped to 0
        assert_eq!(result.uncertainty_lower, Some(0.0));
        assert!(
            (result.uncertainty_upper.unwrap() - (std_expected / 2.0 * 100.0) as f32).abs()
                < TEST_EPSILON_F32
        );
    }

    #[test]
    fn test_calculate_msi_metrics_upper_bound_clamped_at_total() {
        let v1 = make_variant(0.1, 0.8);
        let v2 = make_variant(0.1, 0.8);
        let filtered = FilteredRegions {
            variants: vec![&v1, &v2],
            region_starts: vec![0, 1],
        };
        let output_req = OutputRequirements {
            needs_pseudotime: true,
            needs_distribution: false,
            needs_heatmap: false,
        };

        let result = calculate_msi_metrics(&filtered, 2, 0.5, output_req, 0.05);

        // p_stable = [0.1, 0.1]: P(0) = 0.01, P(1) = 0.18, P(2) = 0.81, so k_map = 2
        // std around k_map = sqrt(2² × 0.01 + 1² × 0.18) = sqrt(0.22)
        let std_expected = (4.0 * 0.01 + 0.18_f64).sqrt();

        assert_eq!(result.k_map, Some(2));
        // k_map + std > 2 regions, so the upper bound is clamped to 2/2 × 100
        assert_eq!(result.uncertainty_upper, Some(100.0));
        assert!(
            (result.uncertainty_lower.unwrap() - ((2.0 - std_expected) / 2.0 * 100.0) as f32).abs()
                < TEST_EPSILON_F32
        );
    }

    /* =========== get_window_slice tests ============ */

    #[test]
    fn test_get_window_slice_basic() {
        let regions = vec![
            make_region("chr1", 100, vec![], true),
            make_region("chr1", 500, vec![], true),
            make_region("chr1", 1100, vec![], true),
            make_region("chr10", 10, vec![], true),
            make_region("chr2", 100, vec![], true),
        ];
        let slice = get_window_slice(&regions, "chr1", 0, 1000);
        assert_eq!(slice.len(), 2);
        assert_eq!(slice[0].start, 100);
        assert_eq!(slice[1].start, 500);

        let slice = get_window_slice(&regions, "chr1", 1000, 2000);
        assert_eq!(slice.len(), 1);
        assert_eq!(slice[0].start, 1100);

        let slice = get_window_slice(&regions, "chr2", 0, 1000);
        assert_eq!(slice.len(), 1);
        assert_eq!(slice[0].start, 100);

        let slice = get_window_slice(&regions, "chr10", 0, 100);
        assert_eq!(slice.len(), 1);
        assert_eq!(slice[0].chrom, "chr10");
        assert_eq!(slice[0].start, 10);
    }

    #[test]
    fn test_get_window_slice_empty_window() {
        let regions = vec![
            make_region("chr1", 100, vec![], true),
            make_region("chr1", 300, vec![], true),
            make_region("chr10", 10, vec![], true),
            make_region("chr2", 100, vec![], true),
        ];

        // Window entirely after all data on the queried chrom
        let slice = get_window_slice(&regions, "chr1", 5000, 6000);
        assert_eq!(slice.len(), 0, "window after all data should be empty");

        // Window entirely before all data on the queried chrom
        let slice = get_window_slice(&regions, "chr1", 0, 50);
        assert_eq!(slice.len(), 0, "window before all data should be empty");

        // Window in a numeric gap between two regions on the same chrom
        let slice = get_window_slice(&regions, "chr1", 150, 250);
        assert_eq!(slice.len(), 0, "window in a numeric gap should be empty");

        // Chromosome missing, but sorting in the middle of the array
        // "chr15" sorts between "chr10" and "chr2" lexicographically.
        let slice = get_window_slice(&regions, "chr15", 0, 1000);
        assert_eq!(
            slice.len(),
            0,
            "chrom missing from the middle of the sort order should be empty"
        );
    }

    #[test]
    fn test_get_window_slice_chrom_not_found() {
        let regions = vec![make_region("chr1", 100, vec![], true)];
        let slice = get_window_slice(&regions, "chr99", 0, 1000);
        assert_eq!(slice.len(), 0);
    }

    #[test]
    fn test_get_window_slice_all_regions_in_window() {
        let regions = vec![
            make_region("chr1", 100, vec![], true),
            make_region("chr1", 200, vec![], true),
        ];
        let slice = get_window_slice(&regions, "chr1", 0, 1_000_000);
        assert_eq!(slice.len(), 2);
    }

    #[test]
    fn test_get_window_slice_boundary_inclusive_exclusive() {
        let regions = vec![
            make_region("chr1", 0, vec![], true),
            make_region("chr1", 1000, vec![], true),
            make_region("chr1", 2000, vec![], true),
        ];

        // window_start is inclusive: region at exactly window_start=1000 included
        let slice = get_window_slice(&regions, "chr1", 1000, 1500);
        assert_eq!(slice.len(), 1);
        assert_eq!(slice[0].start, 1000);

        // window_end is exclusive: region at exactly window_end=1000 excluded
        let slice = get_window_slice(&regions, "chr1", 0, 1000);
        assert_eq!(slice.len(), 1);
        assert_eq!(slice[0].start, 0);
    }

    /* ======== run_windowed_analysis tests ========== */

    #[test]
    fn test_run_windowed_analysis_basic() {
        let regions = vec![
            make_region("chr1", 100, vec![make_variant(0.1, 0.8)], true),
            make_region("chr1", 200, vec![make_variant(0.2, 0.9)], true),
            make_region("chr2", 100, vec![make_variant(0.3, 0.7)], true),
        ];
        let results = run_windowed_analysis(&regions, 0.1, 1_000_000);

        // chr1: 1 window, chr2: 1 window = 2 results
        assert_eq!(results.len(), 2);

        // chr1: 2 regions, p_stable=[0.1, 0.2]
        // P(2)=0.9×0.8=0.72, k_map=2, score=2/2×100=100%, posterior=0.72
        assert_eq!(results[0].chrom, "chr1");
        assert_eq!(results[0].window_start, 0);
        assert_eq!(results[0].regions_in_window, 2);
        assert!((results[0].msi_score - 100.0).abs() < TEST_EPSILON_F32);
        assert!((results[0].posterior_probability - 0.72).abs() < TEST_EPSILON);

        // chr2: 1 region, p_stable=0.3
        // P(1)=0.7 > P(0)=0.3, k_map=1, score=1/1×100=100%, posterior=0.7
        assert_eq!(results[1].chrom, "chr2");
        assert_eq!(results[1].regions_in_window, 1);
        assert!((results[1].msi_score - 100.0).abs() < TEST_EPSILON_F32);
        assert!((results[1].posterior_probability - 0.7).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_windowed_analysis_empty_window_excluded() {
        let regions = vec![
            make_region("chr1", 100, vec![make_variant(0.1, 0.8)], true),
            make_region("chr1", 5_000_100, vec![make_variant(0.2, 0.9)], true),
        ];
        let results = run_windowed_analysis(&regions, 0.1, 1_000_000);

        // Only 2 non-empty windows: 0 and 5_000_000
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].window_start, 0);
        assert_eq!(results[1].window_start, 5_000_000);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "regions=[] guaranteed unreachable")]
    fn test_run_windowed_analysis_empty_regions() {
        let regions: Vec<RegionSummary> = vec![];
        run_windowed_analysis(&regions, 0.1, 1_000_000);
    }

    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "window_size=0 guaranteed unreachable")]
    fn test_run_windowed_analysis_zero_window_size() {
        let regions = vec![make_region("chr1", 100, vec![make_variant(0.1, 0.8)], true)];
        run_windowed_analysis(&regions, 0.1, 0);
    }

    #[test]
    fn test_run_windowed_analysis_af_filter_applied() {
        // All variants have af=0.05, below threshold 0.1 - filtered out
        // Region exists but no variants pass, msi_score=0, posterior=1.0
        let regions = vec![make_region(
            "chr1",
            100,
            vec![make_variant(0.1, 0.05)],
            true,
        )];
        let results = run_windowed_analysis(&regions, 0.1, 1_000_000);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].msi_score, 0.0);
        assert!((results[0].posterior_probability - 1.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_windowed_analysis_uses_per_window_denominator() {
        // 2 regions in window, both unstable, k_map=2, score=2/2×100=100%
        // Windowed uses per-window denominator
        let regions = vec![
            make_region("chr1", 100, vec![make_variant(0.01, 0.8)], true),
            make_region("chr1", 200, vec![make_variant(0.01, 0.9)], true),
        ];
        let results = run_windowed_analysis(&regions, 0.0, 1_000_000);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].regions_in_window, 2);
        assert!((results[0].msi_score - 100.0).abs() < TEST_EPSILON_F32);
    }

    /* ========= run_global_analysis tests =========== */

    #[test]
    fn test_run_global_analysis_basic() {
        let regions = vec![make_region_simple(vec![make_variant(0.1, 0.5)])];
        let output_req = OutputRequirements {
            needs_pseudotime: true,
            needs_distribution: false,
            needs_heatmap: false,
        };
        let results = run_global_analysis(&regions, 100, &[0.0], 0.05, output_req);

        assert_eq!(results.len(), 1);
        let r = &results["0.00"];
        assert_eq!(r.k_map.unwrap(), 1);
        assert!((r.msi_score_map.unwrap() - 1.0).abs() < TEST_EPSILON_F32);
        assert_eq!(r.regions_af_pass.unwrap(), 1);
    }

    /* ====== run_af_evolution_analysis tests ======== */

    #[test]
    fn test_run_af_evolution_analysis_basic() {
        let regions = vec![make_region_simple(vec![make_variant(0.1, 0.5)])];

        let output_req = OutputRequirements {
            needs_pseudotime: false,
            needs_distribution: true,
            needs_heatmap: false,
        };

        let config = AnalysisConfig {
            total_regions: 100,
            sample: "sample1",
            msi_high_threshold: 3.5,
            af_thresholds: vec![0.0],
            num_threads: Some(1),
            window_size: 1_000_000,
            distribution_af: 0.00,
            windowed_af: 0.05,
        };

        let (results, window_results) =
            run_af_evolution_analysis(&regions, config, output_req).unwrap();

        assert!(window_results.is_empty());
        let result = &results["0.00"];

        assert!(result.k_map.is_none());
        assert!(result.msi_score_map.is_none());
        assert!(result.regions_af_pass.is_none());
        assert!(result.uncertainty_lower.is_none());
        assert!(result.uncertainty_upper.is_none());
        assert!(result.map_std_dev.is_none());

        // Distribution SHOULD exist (needs_distribution: true AND af=0.0)
        assert!(result.distribution.is_some());
        let dist = result.distribution.as_ref().unwrap();
        assert_eq!(dist.len(), 2);
        let prob_sum: f64 = dist.iter().map(|d| d.probability).sum();
        assert!((prob_sum - 1.0).abs() < TEST_EPSILON);
        // With prob_absent=0.1, expect P(k=1) > P(k=0)
        assert!(dist[1].probability > dist[0].probability);
    }

    #[test]
    fn test_run_af_evolution_analysis_multiple_af_thresholds() {
        let regions = vec![
            make_region("chr1", 100, vec![make_variant(0.1, 0.9)], true),
            make_region("chr1", 200, vec![make_variant(0.05, 0.5)], true),
        ];

        let output_req = OutputRequirements {
            needs_pseudotime: true,
            needs_distribution: false,
            needs_heatmap: false,
        };

        let config = AnalysisConfig {
            total_regions: 100,
            sample: "sample1",
            msi_high_threshold: 3.5,
            af_thresholds: vec![0.0, 0.5, 1.0],
            num_threads: Some(1),
            window_size: 1_000_000,
            distribution_af: 0.05,
            windowed_af: 0.05,
        };

        let (results, window_results) =
            run_af_evolution_analysis(&regions, config, output_req).unwrap();

        assert!(window_results.is_empty());
        assert_eq!(results.len(), 3, "Should have 3 AF threshold results");

        // AF=1.0: no variants pass (0.9 < 1.0, 0.5 < 1.0)
        let af_1_0 = &results["1.00"];
        assert_eq!(
            af_1_0.regions_af_pass.unwrap(),
            0,
            "No regions should pass AF=1.0 threshold"
        );
        assert_eq!(
            af_1_0.k_map.unwrap(),
            0,
            "k_map should be 0 with no regions"
        );
        assert_eq!(af_1_0.msi_score_map.unwrap(), 0.0, "MSI score should be 0%");
        assert!(
            af_1_0.uncertainty_lower.is_some(),
            "Uncertainty should exist with pseudotime=true"
        );

        // AF=0.5: both variants pass (0.9 ≥ 0.5, 0.5 ≥ 0.5)
        let af_0_5 = &results["0.50"];
        assert_eq!(
            af_0_5.regions_af_pass.unwrap(),
            2,
            "Both regions should pass AF=0.5 threshold"
        );
        assert_eq!(
            af_0_5.k_map.unwrap(),
            2,
            "k_map should be 2 with high instability"
        );
        assert!(
            (af_0_5.msi_score_map.unwrap() - 2.0).abs() < TEST_EPSILON_F32,
            "MSI score should be 2.0%, got {}",
            af_0_5.msi_score_map.unwrap()
        );
        assert!(af_0_5.uncertainty_lower.is_some());
        assert!(af_0_5.uncertainty_upper.is_some());

        // AF=0.0: both variants pass (all pass)
        let af_0_0 = &results["0.00"];
        assert_eq!(
            af_0_0.regions_af_pass.unwrap(),
            2,
            "Both regions should pass AF=0.0 threshold"
        );
        assert_eq!(af_0_0.k_map.unwrap(), 2);
        assert!((af_0_0.msi_score_map.unwrap() - 2.0).abs() < TEST_EPSILON_F32);
        assert!(af_0_0.uncertainty_lower.is_some());

        // Verify no distribution (needs_distribution=false)
        assert!(af_1_0.distribution.is_none());
        assert!(af_0_5.distribution.is_none());
        assert!(af_0_0.distribution.is_none());

        assert_eq!(
            af_0_5.msi_score_map, af_0_0.msi_score_map,
            "MSI scores should match when same regions pass"
        );
    }

    #[test]
    fn test_run_af_evolution_analysis_heatmap_path() {
        let regions = vec![
            make_region("chr1", 100, vec![make_variant(0.1, 0.8)], true),
            make_region("chr1", 200, vec![make_variant(0.2, 0.9)], true),
        ];

        let output_req = OutputRequirements {
            needs_pseudotime: false,
            needs_distribution: false,
            needs_heatmap: true,
        };

        let config = AnalysisConfig {
            total_regions: 100,
            sample: "sample1",
            msi_high_threshold: 3.5,
            af_thresholds: vec![0.0],
            num_threads: Some(1),
            window_size: 1_000_000,
            distribution_af: 0.05,
            windowed_af: 0.05,
        };

        let (global_results, window_results) =
            run_af_evolution_analysis(&regions, config, output_req).unwrap();

        assert!(global_results.is_empty());

        let windows = window_results;
        assert_eq!(windows.len(), 1);
        assert_eq!(windows[0].chrom, "chr1");
        assert_eq!(windows[0].window_start, 0);
        assert_eq!(windows[0].regions_in_window, 2);
        assert!((windows[0].msi_score - 100.0).abs() < TEST_EPSILON_F32);
        assert!((windows[0].posterior_probability - 0.72).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_af_evolution_analysis_nothing_requested() {
        let regions = vec![make_region_simple(vec![make_variant(0.1, 0.5)])];

        let output_req = OutputRequirements {
            needs_pseudotime: false,
            needs_distribution: false,
            needs_heatmap: false,
        };

        let config = AnalysisConfig {
            total_regions: 100,
            sample: "sample1",
            msi_high_threshold: 3.5,
            af_thresholds: vec![0.0],
            num_threads: Some(1),
            window_size: 1_000_000,
            distribution_af: 0.05,
            windowed_af: 0.05,
        };

        let (global_results, window_results) =
            run_af_evolution_analysis(&regions, config, output_req).unwrap();

        assert!(global_results.is_empty());
        assert!(window_results.is_empty());
    }
}

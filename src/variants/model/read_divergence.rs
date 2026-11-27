//! Mapping of reads from their own mismatches.
//!
//! MAPQ only reflects alternative placements within the reference.
//! A read may also come from a source the reference lacks
//! (a paralog, pseudogene or contaminant, sequence that is absent or divergent in the reference).
//! Such a read shows mismatches against the reference at about the divergence of its source,
//! while a read from the mapped locus deviates only by its base errors
//! and by the variation that the other reads of the locus share.
//! For consensus reads, whose base errors are far below any divergence,
//! this makes the read's own mismatches the stronger mapping evidence;
//! for raw reads the base qualities keep the evidence weak,
//! which is why the model applies to consensus reads by default.
//! Preprocessing records the mismatches (see evidence::observations::mismatches); the model runs at calling time,
//! where the scenario tells the divergence and what variation the sample may carry.

use bio::stats::LogProb;
use rgsl::randist::beta::beta_P;
use rgsl::randist::binomial::binomial_pdf;

use crate::grammar;
use crate::variants::evidence::observations::mismatches::Mismatches;
use crate::variants::evidence::observations::pileup::Pileup;

/// Divergence of sources the reference lacks, if the scenario gives none.
pub(crate) const DEFAULT_DIVERGENCE: f64 = 0.03;

/// Reads to which the divergence model applies.
#[derive(
    Display,
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    EnumString,
    EnumIter,
    IntoStaticStr,
    VariantNames,
    Default,
)]
#[strum(serialize_all = "kebab_case")]
#[serde(rename_all = "kebab-case")]
pub enum ReadDivergenceReads {
    /// Consensus reads (see mismatches::is_consensus_read).
    #[default]
    Consensus,
    /// All reads.
    All,
    /// No reads.
    None,
}

/// Settings of the divergence model.
#[derive(Debug, Clone, Copy, Default)]
pub struct ReadDivergenceSettings {
    pub reads: ReadDivergenceReads,
}

/// What the scenario knows about the variation a sample's reads may carry, at a contig.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct SampleVariation {
    /// Ploidy of the sample at the contig, if known.
    pub ploidy: Option<u32>,
    /// Somatic effective mutation rate, if known.
    pub somatic_rate: Option<f64>,
    /// Clustering of somatic mutations, if given.
    pub somatic_clustering: Option<grammar::SomaticClustering>,
}

/// A base of quality 0 is a random base: the error cannot exceed three quarters.
const MAX_ERROR: f64 = 0.75;

fn ln_add_exp(a: f64, b: f64) -> f64 {
    *LogProb(a).ln_add_exp(LogProb(b))
}

/// Rate of mismatches that no other fragment of the locus shares, per aligned base, over the given summaries,
/// by Laplace's rule (one pseudo observation of each kind).
fn background<'a>(summaries: impl Iterator<Item = &'a Mismatches>) -> f64 {
    let (mut bases, mut unshared) = (0u64, 0u64);
    for summary in summaries {
        bases += summary.bases as u64;
        unshared += summary
            .mismatches
            .iter()
            .filter(|(_, shared, _)| *shared == 0)
            .count() as u64;
    }
    (unshared + 1) as f64 / (bases + 2) as f64
}

/// Probability that a base shown by `shared` of `covering` other fragments is established variation of the sample:
/// its count fits a germline VAF of the ploidy (1 / ploidy to 1, equally likely)
/// at least as well as a frequency below the lowest one (flat).
/// The other fragments err at the background rate.
/// Beyond float precision, a frequency far above the lowest germline VAF counts as established.
fn established(shared: u16, covering: u16, ploidy: u32, background: f64) -> f64 {
    let (s, c, p) = (shared as f64, covering as f64, ploidy as f64);
    let germline = (1..=ploidy)
        .map(|j| {
            let vaf = j as f64 / p;
            let shown = vaf * (1.0 - background) + (1.0 - vaf) * background / 3.0;
            binomial_pdf(shared as u32, shown, covering as u32)
        })
        .sum::<f64>()
        / p;
    let below = p * beta_P(1.0 / p, s + 1.0, c - s + 1.0) / (c + 1.0);
    if germline + below > 0.0 {
        germline / (germline + below)
    } else {
        1.0
    }
}

/// Log likelihood ratio of the read's bases under a source of the given divergence against the mapped locus.
///
/// METHOD: both hypotheses read the molecule's true base with the same error:
/// the base quality, or at least the background of unshared mismatches, and at most that of a random base.
/// An error turns the base into each of the three other bases with a third of that rate.
/// So a mismatch of low quality says nothing either way.
/// Under another source, every non-reference base of the molecule comes from the divergence.
/// Under the mapped locus, a base that no other fragment shares is an error.
/// A base that other fragments share is established variation of the sample (see `established`),
/// which the molecule carries with the frequency of the other fragments, or an error.
/// Rarer shared variation is somatic and has to exist first: as independent mutations it is rarer than errors,
/// but as a cluster it is not (see grammar::SomaticClustering):
/// a cluster covers the molecule with probability fraction * somatic rate / cluster rate
/// and then mutates each base at the cluster rate, so the shared bases that are not established come from it.
/// The positions of the mismatches are not recorded, so the matching bases are taken as outside the cluster.
/// Bases that no other fragment shares stay errors even then: a cluster in a single molecule looks like a divergent read.
/// Each read is scored on its own, so a group of reads from one subclone pays for the subclone once per read;
/// and a dense cluster of low frequency looks like a paralog of low divergence, which the cluster rate against the divergence decides.
/// Without a ploidy, every shared base counts as established.
/// Mismatches beyond the stored ones count as uninformative.
fn log_likelihood_ratio(
    summary: &Mismatches,
    background: f64,
    divergence: f64,
    variation: &SampleVariation,
) -> f64 {
    let d = divergence;
    let prob_other = |error: f64| {
        d / 3.0 * (1.0 - error) + (1.0 - d) * error / 3.0 + 2.0 * d / 3.0 * error / 3.0
    };
    let cluster = match (variation.somatic_clustering, variation.somatic_rate) {
        (Some(clustering), Some(rate)) => Some((
            clustering.fraction * rate / clustering.rate,
            clustering.rate,
        )),
        _ => None,
    };
    let mut log_other = 0.0;
    // the shared bases as errors, and as a cluster
    let mut log_mapped = 0.0;
    let mut log_clustered = 0.0;
    for (qual, shared, covering) in &summary.mismatches {
        let error = 10f64
            .powf(-(*qual as f64) / 10.0)
            .max(background)
            .min(MAX_ERROR);
        log_other += prob_other(error).ln();
        let (established, frequency) = if *shared == 0 || *covering == 0 {
            (0.0, 0.0)
        } else {
            let frequency = *shared as f64 / *covering as f64;
            match variation.ploidy {
                Some(ploidy) => (
                    established(*shared, *covering, ploidy, background),
                    frequency,
                ),
                None => (1.0, frequency),
            }
        };
        log_mapped += (established * frequency + (1.0 - established) * error / 3.0).ln();
        if let Some((_, rate)) = cluster {
            let mutated = if *shared == 0 { error } else { rate };
            log_clustered += (established * frequency + (1.0 - established) * mutated / 3.0).ln();
        }
    }
    let matches = summary.bases.saturating_sub(summary.len() as u32) as f64;
    log_other += matches * (1.0 - d).ln();
    if let Some((covered, _)) = cluster {
        log_mapped = ln_add_exp(log_mapped, covered.ln() + log_clustered);
    }
    log_other - log_mapped
}

/// Update the mapping probability of each observation the model applies to with the evidence of its mismatches,
/// and let a mismapped read support the alt allele only where its source differs from the mapped locus.
pub(crate) fn adjust_read_divergence(
    pileup: &mut Pileup,
    divergence: f64,
    settings: &ReadDivergenceSettings,
    variation: &SampleVariation,
) {
    let applies = |consensus: bool| match settings.reads {
        ReadDivergenceReads::Consensus => consensus,
        ReadDivergenceReads::All => true,
        ReadDivergenceReads::None => false,
    };
    let background = background(
        pileup
            .read_observations()
            .iter()
            .filter(|obs| applies(obs.consensus) && obs.mismatches.bases > 0)
            .map(|obs| &obs.mismatches),
    );
    for obs in pileup.read_observations_mut().iter_mut() {
        if !applies(obs.consensus) || obs.mismatches.bases == 0 {
            continue;
        }
        let log_ratio = log_likelihood_ratio(&obs.mismatches, background, divergence, variation);
        let prob_mapping = obs.prob_mapping();
        let mismapped = LogProb(*prob_mapping.ln_one_minus_exp() + log_ratio);
        let prob_mismapping = LogProb(*mismapped - *mismapped.ln_add_exp(prob_mapping));
        // METHOD: another source shows the alt base where it differs from the mapped locus.
        let prob_missed_allele = if obs.prob_alt() > obs.prob_ref() {
            LogProb((divergence / 3.0).ln())
        } else {
            LogProb((1.0 - divergence).ln())
        };
        obs.set_mapping(prob_mismapping.ln_one_minus_exp(), prob_missed_allele);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn summary(bases: u32, mismatches: Vec<(u8, u16, u16)>) -> Mismatches {
        Mismatches {
            bases,
            mismatches,
            more: 0,
        }
    }

    // background of consensus reads
    const BG: f64 = 4e-5;
    const D: f64 = 0.03;

    fn ratio(summary: &Mismatches, variation: &SampleVariation) -> f64 {
        log_likelihood_ratio(summary, BG, D, variation)
    }

    const NONE: SampleVariation = SampleVariation {
        ploidy: None,
        somatic_rate: None,
        somatic_clustering: None,
    };
    const DIPLOID: SampleVariation = SampleVariation {
        ploidy: Some(2),
        somatic_rate: Some(1e-6),
        somatic_clustering: None,
    };
    const CLUSTERED: SampleVariation = SampleVariation {
        somatic_clustering: Some(grammar::SomaticClustering {
            fraction: 0.05,
            rate: 0.02,
        }),
        ..DIPLOID
    };

    #[test]
    fn test_unshared_mismatches() {
        // a clean read is less likely from another source
        let clean = ratio(&summary(60, vec![]), &NONE);
        assert!(clean < 0.0);
        // unexplained mismatches of high quality point to another source, more so than the clean bases speak against it
        let high = ratio(&summary(60, vec![(90, 0, 20); 3]), &NONE);
        assert!(high > 2f64.ln(), "{}", high);
        // mismatches of low quality say little either way: a tenth of the high quality ones at most
        let low = ratio(&summary(60, vec![(10, 0, 20); 3]), &NONE);
        assert!(
            low - clean < 0.1 * (high - clean),
            "{} {} {}",
            low,
            clean,
            high
        );
        // a base of quality 0 is finite and uninformative: the read counts as one base shorter
        let zero = ratio(&summary(60, vec![(0, 0, 20)]), &NONE);
        let shorter = ratio(&summary(59, vec![]), &NONE);
        assert!(
            zero.is_finite() && (zero - shorter).abs() < 1e-6,
            "{} {}",
            zero,
            shorter
        );
        // unshared mismatches stay errors with clustering: a single divergent read is still seen
        let divergent = summary(150, vec![(90, 0, 170); 3]);
        assert!((ratio(&divergent, &CLUSTERED) - ratio(&divergent, &DIPLOID)).abs() < 1e-4);
    }

    #[test]
    fn test_established() {
        // half of the other fragments or nearly all of them: germline of a diploid
        assert!(established(90, 180, 2, BG) > 0.8);
        assert!(established(178, 180, 2, BG) > 0.9);
        // a few of them: not
        assert!(established(1, 173, 2, BG) < 1e-3);
        assert!(established(60, 180, 2, BG) < 0.01);
        // a third of the fragments is germline for a triploid only
        assert!(established(60, 180, 3, BG) > 0.5);
        // a count somewhat below half still passes mostly, and the weight falls off smoothly
        assert!(established(76, 170, 2, BG) > 0.4);
        assert!(established(72, 170, 2, BG) < established(76, 170, 2, BG));
        assert!(established(76, 170, 2, BG) < established(85, 170, 2, BG));
        // the degenerate frequencies at float precision
        assert!(established(7500, 10000, 2, BG) > 0.9);
        assert!(established(4500, 10000, 2, BG) < 0.01);
    }

    #[test]
    fn test_shared_mismatches() {
        let clean = ratio(&summary(96, vec![]), &NONE);
        // bases shared by half the fragments are established and speak for the mapped locus
        let het = summary(96, vec![(90, 86, 173), (90, 88, 173), (90, 85, 173)]);
        assert!(
            ratio(&het, &DIPLOID) < clean,
            "{} {}",
            ratio(&het, &DIPLOID),
            clean
        );
        // without a ploidy every shared base is established, which makes little difference here
        assert!(ratio(&het, &NONE) < clean);
        assert!((ratio(&het, &NONE) - ratio(&het, &DIPLOID)).abs() < 1.0);
        // one base shared by one other fragment is an error, with or without clustering
        let single = summary(96, vec![(90, 1, 173)]);
        assert!((ratio(&single, &CLUSTERED) - ratio(&single, &DIPLOID)).abs() < 1e-2);
        // three bases shared by one other fragment: another source, unless somatic mutations cluster
        let triplet = summary(96, vec![(90, 1, 173); 3]);
        let independent = ratio(&triplet, &DIPLOID);
        let in_cluster = ratio(&triplet, &CLUSTERED);
        assert!(independent > 1000f64.ln(), "{}", independent);
        assert!(
            in_cluster < independent - 10f64.ln(),
            "{} {}",
            in_cluster,
            independent
        );
        // five bases shared by four other fragments: another source, and less certain so with clustering
        let group = summary(96, vec![(90, 4, 300); 5]);
        assert!(ratio(&group, &DIPLOID) > independent);
        assert!(ratio(&group, &CLUSTERED) < ratio(&group, &DIPLOID));
        assert!(
            ratio(&group, &CLUSTERED) > 100f64.ln(),
            "{}",
            ratio(&group, &CLUSTERED)
        );
    }

    #[test]
    fn test_background() {
        let summaries = [
            summary(100, vec![(90, 0, 5), (90, 3, 5)]),
            summary(100, vec![]),
        ];
        assert!((background(summaries.iter()) - 2.0 / 202.0).abs() < 1e-12);
    }
}

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
//! where the divergence comes from the scenario.

use bio::stats::LogProb;

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
#[derive(Debug, Clone, Copy)]
pub struct ReadDivergenceSettings {
    pub reads: ReadDivergenceReads,

    /// A base shared by at least this fraction of the other fragments covering the position is variation of the sample.
    /// A rarer shared base is divergence, since a few reads from another source share their mismatches as well.
    /// Set below one over the ploidy of the sample; the default suits diploid samples with some allelic imbalance.
    /// This is a stand-in, see log_likelihood_ratio().
    pub variation_min_frequency: f64,
}

impl Default for ReadDivergenceSettings {
    fn default() -> Self {
        ReadDivergenceSettings {
            reads: ReadDivergenceReads::Consensus,
            variation_min_frequency: 0.2,
        }
    }
}

/// A base of quality 0 is a random base: the error cannot exceed three quarters.
const MAX_ERROR: f64 = 0.75;

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

/// Log likelihood ratio of the read's bases under a source of the given divergence against the mapped locus.
///
/// METHOD: both hypotheses read the molecule's true base with the same error:
/// the base quality, or at least the background of unshared mismatches, and at most that of a random base.
/// An error turns the base into each of the three other bases with a third of that rate.
/// The hypotheses differ in the true base: under the mapped locus it is a non-reference base
/// with the frequency of the sample's variation (see variation_min_frequency), under another source with the divergence.
/// So a mismatch of low quality says nothing either way.
/// Matching bases count the divergence against the sample's variation, taken as absent there;
/// mismatches beyond the stored ones count as uninformative.
/// TODO: replace the minimum frequency by a mixture of haplotypes per locus, those of the sample against those of
/// another source, with the scenario's heterozygosity, somatic mutation rate and ploidy and the read length as priors
/// (a long read may well carry several somatic differences).
fn log_likelihood_ratio(
    summary: &Mismatches,
    background: f64,
    divergence: f64,
    variation_min_frequency: f64,
) -> f64 {
    let mut log_ratio = 0.0;
    for (qual, shared, covering) in &summary.mismatches {
        let error = 10f64
            .powf(-(*qual as f64) / 10.0)
            .max(background)
            .min(MAX_ERROR);
        let frequency = if *covering > 0 {
            *shared as f64 / *covering as f64
        } else {
            0.0
        };
        let variation = if frequency >= variation_min_frequency {
            frequency
        } else {
            0.0
        };
        let prob_mapped = variation * (1.0 - error) + (1.0 - variation) * error / 3.0;
        let prob_other = divergence / 3.0 * (1.0 - error)
            + (1.0 - divergence) * error / 3.0
            + 2.0 * divergence / 3.0 * error / 3.0;
        log_ratio += prob_other.ln() - prob_mapped.ln();
    }
    let matches = summary.bases.saturating_sub(summary.len() as u32);
    log_ratio + matches as f64 * (1.0 - divergence).ln()
}

/// Update the mapping probability of each observation the model applies to with the evidence of its mismatches,
/// and let a mismapped read support the alt allele only where its source differs from the mapped locus.
pub(crate) fn adjust_read_divergence(
    pileup: &mut Pileup,
    divergence: f64,
    settings: &ReadDivergenceSettings,
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
        let log_ratio = log_likelihood_ratio(
            &obs.mismatches,
            background,
            divergence,
            settings.variation_min_frequency,
        );
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

    #[test]
    fn test_log_likelihood_ratio() {
        let (bg, d, floor) = (1e-3, 0.03, 0.2);
        // a clean read is less likely from another source
        let clean = log_likelihood_ratio(&summary(60, vec![]), bg, d, floor);
        assert!(clean < 0.0);
        // unexplained mismatches of high quality point to another source, more so than the clean bases speak against it
        let high = log_likelihood_ratio(&summary(60, vec![(90, 0, 20); 3]), bg, d, floor);
        assert!(high > 2f64.ln(), "{}", high);
        // mismatches of low quality say little either way: a tenth of the high quality ones at most
        let low = log_likelihood_ratio(&summary(60, vec![(10, 0, 20); 3]), bg, d, floor);
        assert!(
            low - clean < 0.1 * (high - clean),
            "{} {} {}",
            low,
            clean,
            high
        );
        // a base of quality 0 is finite and uninformative: the read counts as one base shorter
        let zero = log_likelihood_ratio(&summary(60, vec![(0, 0, 20)]), bg, d, floor);
        let shorter = log_likelihood_ratio(&summary(59, vec![]), bg, d, floor);
        assert!(
            zero.is_finite() && (zero - shorter).abs() < 1e-6,
            "{} {}",
            zero,
            shorter
        );
        // a mismatch shared by most other fragments is variation of the sample and speaks for the mapped locus
        let shared = log_likelihood_ratio(&summary(60, vec![(90, 10, 20)]), bg, d, floor);
        assert!(shared < clean, "{} {}", shared, clean);
        // a mismatch shared by a few fragments counts as divergence
        let rare = log_likelihood_ratio(&summary(60, vec![(90, 1, 20)]), bg, d, floor);
        assert!(rare > clean, "{} {}", rare, clean);
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

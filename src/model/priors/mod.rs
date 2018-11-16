pub mod infinite_sites_neutral_variation;
pub mod single_cell_bulk;
pub mod tumor_normal;
// TODO disable Tumor-Normal-Relapse model until i
//pub mod tumor_normal_relapse;
pub mod flat;
pub mod normal;

use bio::stats::LogProb;

use model::{AlleleFreq, AlleleFreqs, Variant};

pub use priors::infinite_sites_neutral_variation::InfiniteSitesNeutralVariationModel;
pub use priors::single_cell_bulk::SingleCellBulkModel;
pub use priors::tumor_normal::TumorNormalModel;
//pub use priors::tumor_normal_relapse::TumorNormalRelapseModel;
pub use model::PairPileup;
pub use priors::flat::FlatNormalNormalModel;
pub use priors::flat::FlatTumorNormalModel;

/// A prior model of the allele frequency spectrum.
pub trait Model<A: AlleleFreqs> {
    /// Calculate prior probability of given allele frequency.
    fn prior_prob(&self, af: AlleleFreq, variant: &Variant) -> LogProb;

    /// Return allele frequency spectrum.
    fn allele_freqs(&self) -> &A;

    fn marginal_prob<L>(&self, likelihood: &L, variant: &Variant, n_obs: usize) -> LogProb
    where
        L: Fn(AlleleFreq) -> LogProb;

    fn joint_prob<L>(&self, af: &A, likelihood: &L, variant: &Variant, n_obs: usize) -> LogProb
    where
        L: Fn(AlleleFreq) -> LogProb;

    /// Calculate maximum a posteriori probability estimate of allele frequency.
    fn map<L>(&self, likelihood: &L, variant: &Variant, n_obs: usize) -> AlleleFreq
    where
        L: Fn(AlleleFreq) -> LogProb;
}

/// A prior model for sample pairs.
pub trait PairModel<A: AlleleFreqs, B: AlleleFreqs> {
    /// Calculate prior probability of given combination of allele frequencies.
    fn prior_prob(&self, af1: AlleleFreq, af2: AlleleFreq, variant: &Variant) -> LogProb;

    /// Calculate joint probability of prior with likelihoods for given allele frequency ranges.
    fn joint_prob(&self, af1: &A, af2: &B, pileup: &mut PairPileup<A, B, Self>) -> LogProb;

    /// Calculate marginal probability.
    fn marginal_prob(&self, pileup: &mut PairPileup<A, B, Self>) -> LogProb;

    /// Calculate maximum a posteriori probability estimate of allele frequencies.
    fn map(&self, pileup: &mut PairPileup<A, B, Self>) -> (AlleleFreq, AlleleFreq);

    /// Return allele frequency spectra.
    fn allele_freqs(&self) -> (&A, &B);
}

/// A prior model for trios.
pub trait TrioModel<A: AlleleFreqs, B: AlleleFreqs, C: AlleleFreqs> {
    /// Calculate prior probability of given combination of allele frequencies.
    fn prior_prob(
        &self,
        af1: AlleleFreq,
        af2: AlleleFreq,
        af3: AlleleFreq,
        variant: &Variant,
    ) -> LogProb;

    /// Calculate joint probability of prior with likelihoods for given allele frequency ranges.
    fn joint_prob<L, O, Q>(
        &self,
        af1: &A,
        af2: &B,
        af3: &C,
        likelihood1: &L,
        likelihood2: &O,
        likelihood3: &Q,
        variant: &Variant,
        n_obs1: usize,
        n_obs2: usize,
        n_obs3: usize,
    ) -> LogProb
    where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        Q: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb;

    /// Calculate marginal probability.
    fn marginal_prob<L, O, Q>(
        &self,
        likelihood1: &L,
        likelihood2: &O,
        likelihood3: &Q,
        variant: &Variant,
        n_obs1: usize,
        n_obs2: usize,
        n_obs3: usize,
    ) -> LogProb
    where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        Q: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb;

    /// Calculate maximum a posteriori probability estimate of allele frequencies.
    fn map<L, O, Q>(
        &self,
        likelihood1: &L,
        likelihood2: &O,
        likelihood3: &Q,
        variant: &Variant,
        n_obs1: usize,
        n_obs2: usize,
        n_obs3: usize,
    ) -> (AlleleFreq, AlleleFreq, AlleleFreq)
    where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        Q: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb;

    /// Return allele frequency spectra.
    fn allele_freqs(&self) -> (&A, &B, &C);
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::stats::Prob;
    use itertools_num::linspace;
    use model::evidence::observation::Observation;
    use model::tests::observation;
    use model::AlleleFreq;
    use model::Variant;

    pub fn create_obs_vector(n_obs_ref: usize, n_obs_alt: usize) -> Vec<Observation> {
        let obs_ref_abs = observation(LogProb::ln_one(), LogProb::ln_zero(), LogProb::ln_one());
        let obs_alt_abs = observation(LogProb::ln_one(), LogProb::ln_one(), LogProb::ln_zero());

        let mut obs = Vec::new();
        for _ in 0..n_obs_ref {
            obs.push(obs_ref_abs.clone());
        }
        for _ in 0..n_obs_alt {
            obs.push(obs_alt_abs.clone());
        }
        obs
    }

    #[test]
    fn test_infinite_sites_neutral_variation() {
        let ploidy = 2;
        let het = Prob(0.001);
        let model = InfiniteSitesNeutralVariationModel::new(1, ploidy, het);
        assert_relative_eq!(model.prior_prob(1).exp(), 0.001);
        assert_relative_eq!(model.prior_prob(2).exp(), 0.0005);
        assert_relative_eq!(model.prior_prob(0).exp(), 0.9985, epsilon = 0.00001);
    }

    #[test]
    fn test_tumor() {
        let variant = Variant::Deletion(3);
        let model = TumorNormalModel::new(2, 300.0, 1.0, 1.0, 3e9 as u64, Prob(0.001));
        println!(
            "af=0.0,0.0 -> {}",
            *model.prior_prob(AlleleFreq(0.0), AlleleFreq(0.0), &variant)
        );
        println!(
            "af=0.5,0.5 -> {}",
            *model.prior_prob(AlleleFreq(0.5), AlleleFreq(0.5), &variant)
        );

        for af in linspace(0.0, 1.0, 20) {
            println!(
                "normal={} af={} p={}",
                0.0,
                af,
                *model.prior_prob(AlleleFreq(af), AlleleFreq(0.0), &variant)
            );
            println!(
                "normal={} af={} p={}",
                0.5,
                af,
                *model.prior_prob(AlleleFreq(af), AlleleFreq(0.5), &variant)
            );
            println!(
                "normal={} af={} p={}",
                1.0,
                af,
                *model.prior_prob(AlleleFreq(af), AlleleFreq(1.0), &variant)
            );
        }
    }
}

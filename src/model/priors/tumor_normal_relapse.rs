
use bio::stats::{LogProb, Prob};
use itertools::Itertools
use itertools_num::linspace;
use ordered_float::NotNaN;

use priors::{PairModel, TumorNormalModel, InfiniteSitesNeutralVariationModel};
use priors::TrioModel;
use model::{Variant, ContinuousAlleleFreqs, DiscreteAlleleFreqs, AlleleFreq};


pub struct TumorNormalRelapseModel {
    primary_model: TumorNormalModel,
    relapse_model: TumorNormalModel
}


impl TumorNormalRelapseModel {
    /// Create new model.
    ///
    /// # Arguments
    ///
    /// * `ploidy` - the ploidy in the corresponding normal sample (e.g. 2 for diploid)
    /// * `effective_mutation_rate` - the SNV mutation rate per effective cell division in the tumor
    /// * `deletion_factor` - ratio of deletions compared to SNV mutation rate
    /// * `insertion_factor` - ratio of insertions compared to SNV mutation rate
    /// * `genome_size` - the size of the genome
    /// * `heterozygosity` - expected heterozygosity in the corresponding normal
    pub fn new(
        ploidy: u32,
        effective_mutation_rate: f64,
        deletion_factor: f64,
        insertion_factor: f64,
        genome_size: u64,
        heterozygosity: Prob) -> Self {

        TumorNormalRelapseModel {
            primary_model: TumorNormalModel::new(
                ploidy,
                effective_mutation_rate,
                deletion_factor,
                insertion_factor,
                genome_size,
                heterozygosity
            ),
            relapse_model: TumorNormalModel::new(
                ploidy,
                effective_mutation_rate,
                deletion_factor,
                insertion_factor,
                genome_size,
                // this heterozygosity is for the case that we have a variant in the surviving subclone
                // hence, we assume 0 = 1.0 - het * sum_i 1/i and calculate het below
                Prob(1.0 / InfiniteSitesNeutralVariationModel::allele_freq_sum(ploidy))
            )
        }
    }
}


impl TrioModel<ContinuousAlleleFreqs, ContinuousAlleleFreqs, DiscreteAlleleFreqs> for TumorNormalRelapseModel {
    fn prior_prob(&self, af_tumor: AlleleFreq, af_relapse: AlleleFreq, af_normal: AlleleFreq, variant: Variant) -> LogProb {
        let p_tumor_normal = self.primary_model.prior_prob(af_tumor, af_normal, variant);
        // af_tumor = 1 or 0
        let p_relapse = if *af_tumor == 1.0 || *af_tumor == 0.0 {
            let af_relapse_somatic = af_relapse - af_tumor;
            self.relapse_model.somatic_prior_prob(af_relapse_somatic, variant)
        } else {
            let surviving_subclone_afs = self.primary_model.allele_freqs().1;
            // TODO with prob af_tumor I have either 0.5 or 1.0 in the surviving subclone
            let mut summands = surviving_subclone_afs.iter().skip(1).map(|af| self.relapse_model.prior_prob(af_relapse, *af, variant)).collect_vec();
            // TODO with probability 1 - af_tumor, the surviving subclone has allele frequency 0.0
            summands.push(LogProb(af_tumor.ln()).ln_one_minus_exp() + self.relapse_model.somatic_prior_prob(af_relapse, variant));
            LogProb::ln_sum_exp(&summands)
        };

        p_tumor_normal + p_relapse
    }

    fn joint_prob<L, O, Q>(
        &self,
        af_tumor: &ContinuousAlleleFreqs,
        af_relapse: &ContinuousAlleleFreqs,
        af_normal: &DiscreteAlleleFreqs,
        likelihood_tumor: &L,
        likelihood_relapse: &O,
        likelihood_normal: &Q,
        variant: Variant
    ) -> LogProb where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        Q: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
        let prob = LogProb::ln_sum_exp(&af_normal.iter().map(|&af_normal| {
            let tumor_density = |af_tumor| {
                let af_tumor = AlleleFreq(af_tumor);

                let relapse_density = |af_relapse| {
                    let af_relapse = AlleleFreq(af_relapse);
                    self.prior_prob(af_tumor, af_relapse, af_normal, variant) +
                    likelihood_relapse(af_relapse, af_normal)
                };

                LogProb::ln_simpsons_integrate_exp(&relapse_density, *af_relapse.start, *af_relapse.end, self.relapse_model.grid_points) +
                likelihood_tumor(af_tumor, af_normal)
            };

            let p_tumor = if af_tumor.start == af_tumor.end {
                tumor_density(*af_tumor.start)
            } else {
                LogProb::ln_simpsons_integrate_exp(&tumor_density, *af_tumor.start, *af_tumor.end, self.primary_model.grid_points)
            };
            let p_normal = likelihood_normal(af_normal, AlleleFreq(0.0));
            let prob = p_tumor + p_normal;

            prob
        }).collect_vec());

        prob
    }

    fn marginal_prob<L, O, Q>(
        &self,
        likelihood_tumor: &L,
        likelihood_relapse: &O,
        likelihood_normal: &Q,
        variant: Variant
    ) -> LogProb where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        Q: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
        self.joint_prob(
            self.allele_freqs().0,
            self.allele_freqs().1,
            self.allele_freqs().2,
            likelihood_tumor,
            likelihood_relapse,
            likelihood_normal,
            variant
        ).ln_add_exp(
            // add prob for allele frequency zero (the density is non-continuous there)
            // TODO check for further non-continuous spots...
            self.joint_prob(
                &(AlleleFreq(0.0)..AlleleFreq(0.0)),
                &(AlleleFreq(0.0)..AlleleFreq(0.0)),
                &vec![AlleleFreq(0.0)],
                likelihood_tumor,
                likelihood_relapse,
                likelihood_normal,
                variant
            )
        )
    }

    fn map<L, O, Q>(
        &self,
        likelihood_tumor: &L,
        likelihood_relapse: &O,
        likelihood_normal: &Q,
        variant: Variant
    ) -> (AlleleFreq, AlleleFreq, AlleleFreq) where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        Q: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
        let af_tumor = linspace(*self.primary_model.allele_freqs_tumor.start, *self.primary_model.allele_freqs_tumor.end, self.primary_model.grid_points);

        let (_, ((map_normal, map_tumor), map_relapse)) = self.primary_model.allele_freqs().1.iter().cartesian_product(af_tumor.clone()).cartesian_product(af_tumor.clone()).minmax_by_key(
            |&((&af_normal, af_tumor), af_relapse)| {
                let af_tumor = AlleleFreq(af_tumor);
                let af_relapse = AlleleFreq(af_relapse);
                let p = self.prior_prob(af_tumor, af_relapse, af_normal, variant) +
                        likelihood_tumor(af_tumor, af_normal) +
                        likelihood_relapse(af_relapse, af_normal) +
                        likelihood_normal(af_normal, AlleleFreq(0.0));
                //println!("af {} vs {} = {} (prior={} tumor={} normal={})", *af_tumor, af_normal, *p, *self.prior_prob(af_tumor, af_normal, variant), *likelihood_tumor(af_tumor, af_normal), *likelihood_normal(af_normal, AlleleFreq(0.0)));
                NotNaN::new(*p).expect("posterior probability is NaN")
            }
        ).into_option().expect("prior has empty allele frequency spectrum");

        (AlleleFreq(map_tumor), AlleleFreq(map_relapse), *map_normal)
    }

    fn allele_freqs(&self) -> (&ContinuousAlleleFreqs, &ContinuousAlleleFreqs, &DiscreteAlleleFreqs) {
        (self.primary_model.allele_freqs().0, self.relapse_model.allele_freqs().0, self.primary_model.allele_freqs().1)
    }
}

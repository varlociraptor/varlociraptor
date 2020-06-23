// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::fmt::Debug;
use std::rc::Rc;
use std::str;
use std::sync::Arc;

use anyhow::Result;
use bio::stats::{self, pairhmm::PairHMM, LogProb, Prob};
use bio_types::genome::{self, AbstractInterval};
use rust_htslib::bam;

use crate::variants::evidence::realignment::edit_distance::EditDistanceCalculation;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, ReferenceEmissionParams};
use crate::variants::types::{AlleleSupport, AlleleSupportBuilder};

pub(crate) mod edit_distance;
pub(crate) mod pairhmm;

pub(crate) trait Realignable<'a> {
    type EmissionParams: stats::pairhmm::EmissionParameters + pairhmm::RefBaseEmission;

    fn alt_emission_params(
        &self,
        read_emission_params: Rc<ReadEmission<'a>>,
        ref_seq: Arc<Vec<u8>>,
        ref_window: usize,
    ) -> Self::EmissionParams;
}

#[derive(Debug, Clone)]
pub(crate) struct Realigner {
    pairhmm: PairHMM,
    gap_params: pairhmm::GapParams,
    max_window: u64,
    ref_seq: Arc<Vec<u8>>,
}

impl Realigner {
    /// Create a new instance.
    pub(crate) fn new(ref_seq: Arc<Vec<u8>>, gap_params: pairhmm::GapParams, max_window: u64) -> Self
where {
        Realigner {
            gap_params,
            pairhmm: PairHMM::new(),
            max_window,
            ref_seq,
        }
    }

    pub(crate) fn allele_support<'a, V>(
        &mut self,
        record: &'a bam::Record,
        locus: &genome::Interval,
        variant: &V,
    ) -> Result<AlleleSupport>
    where
        V: Realignable<'a>,
    {
        let read_seq: bam::record::Seq<'a> = record.seq();
        let read_qual = record.qual();
        let cigar = record.cigar_cached().unwrap();
        let locus_start = locus.range().start;
        let locus_end = locus.range().end;

        let (read_offset, read_end, breakpoint, overlap) = match (
            cigar.read_pos(locus_start as u32, true, true)?,
            cigar.read_pos(locus_end as u32, true, true)?,
        ) {
            // read encloses variant
            (Some(qstart), Some(qend)) => {
                let qstart = qstart as usize;
                // exclusive end of variant
                let qend = qend as usize;
                // ensure that distance between qstart and qend does not make the window too
                // large
                let max_window = (self.max_window as usize).saturating_sub((qend - qstart) / 2);
                let mut read_offset = qstart.saturating_sub(max_window);
                let mut read_end = cmp::min(qend + max_window as usize, read_seq.len());

                // correct for reads that enclose the entire variant while that exceeds the maximum pattern len
                let exceed = (read_end - read_offset)
                    .saturating_sub(EditDistanceCalculation::max_pattern_len());
                if exceed > 0 {
                    read_offset += exceed / 2;
                    read_end -= (exceed as f64 / 2.0).ceil() as usize;
                }
                (read_offset, read_end, locus_start as usize, true)
            }

            (Some(qstart), None) => {
                let qstart = qstart as usize;
                let read_offset = qstart.saturating_sub(self.max_window as usize);
                let read_end = cmp::min(qstart + self.max_window as usize, read_seq.len());
                (read_offset, read_end, locus_start as usize, true)
            }
            (None, Some(qend)) => {
                let qend = qend as usize;
                let read_offset = qend.saturating_sub(self.max_window as usize);
                let read_end = cmp::min(qend + self.max_window as usize, read_seq.len());
                (read_offset, read_end, locus_end as usize, true)
            }
            (None, None) => {
                let m = read_seq.len() / 2;
                let read_offset = m.saturating_sub(self.max_window as usize);
                let read_end = cmp::min(m + self.max_window as usize - 1, read_seq.len());
                let breakpoint = record.pos() as usize + m;
                // The following should only happen with deletions.
                // It occurs if the read comes from ref allele and is mapped within start
                // and end of deletion. Usually, such reads strongly support the ref allele.
                let read_enclosed_by_variant =
                    record.pos() >= locus_start as i64 && cigar.end_pos() <= locus_end as i64;
                (read_offset, read_end, breakpoint, read_enclosed_by_variant)
            }
        };

        // the window on the reference should be a bit larger to allow some flexibility with close
        // indels. But it should not be so large that the read can align outside of the breakpoint.
        let ref_window = (self.max_window as f64 * 1.5) as usize;

        // read emission
        let read_emission = Rc::new(ReadEmission::new(
            read_seq,
            read_qual,
            read_offset,
            read_end,
        ));
        let edit_dist = EditDistanceCalculation::new((read_offset..read_end).map(|i| read_seq[i]));

        if !overlap {
            // If there is no overlap, normalization below would anyway lead to 0.5 vs 0.5,
            // multiplied with certainty estimate. Hence, we can skip the entire HMM calculation!
            let p = LogProb(0.5f64.ln());
            return Ok(AlleleSupportBuilder::default()
                .prob_ref_allele(p)
                .prob_alt_allele(p)
                .no_strand_info()
                .build()
                .unwrap());
        }

        // ref allele
        let mut prob_ref = self.prob_allele(
            ReferenceEmissionParams {
                ref_seq: Arc::clone(&self.ref_seq),
                ref_offset: breakpoint.saturating_sub(ref_window),
                ref_end: cmp::min(breakpoint + ref_window, self.ref_seq.len()),
                read_emission: Rc::clone(&read_emission),
            },
            &edit_dist,
        );

        let mut prob_alt = self.prob_allele(
            variant.alt_emission_params(
                Rc::clone(&read_emission),
                Arc::clone(&self.ref_seq),
                ref_window,
            ),
            &edit_dist,
        );

        assert!(!prob_ref.is_nan());
        assert!(!prob_alt.is_nan());

        // METHOD: Normalize probabilities. By this, we avoid biases due to proximal variants that are in
        // cis with the considered one. They are normalized away since they affect both ref and alt.
        // In a sense, this assumes that the two considered alleles are the only possible ones.
        // However, if the read actually comes from a third allele, both probabilities will be
        // equally bad, and the normalized one will not prefer any of them.
        // This is ok, because for the likelihood function only the ratio between the two
        // probabilities is relevant!

        if prob_ref != LogProb::ln_zero() && prob_alt != LogProb::ln_zero() {
            // METHOD: Only perform normalization if both probs are non-zero
            // otherwise, we would artificially magnify the ratio
            // (compared to an epsilon for the zero case).
            let prob_total = prob_alt.ln_add_exp(prob_ref);
            prob_ref -= prob_total;
            prob_alt -= prob_total;
        }

        if prob_ref == LogProb::ln_zero() && prob_alt == LogProb::ln_zero() {
            // METHOD: if both are zero, use 0.5 instead. Since only the ratio matters, this
            // has the same effect, without rendering the entire pileup likelihood zero.
            prob_ref = LogProb::from(Prob(0.5));
            prob_alt = prob_ref;
            debug!(
                "Record {} neither supports reference nor alternative allele during realignment.",
                str::from_utf8(record.qname()).unwrap()
            );
        }

        let mut builder = AlleleSupportBuilder::default();

        builder.prob_ref_allele(prob_ref).prob_alt_allele(prob_alt);

        if prob_ref != prob_alt {
            builder.register_record(record);
        } else {
            // METHOD: if record is not informative, we don't want to
            // retain its information (e.g. strand).
            builder.no_strand_info();
        }

        Ok(builder.build().unwrap())
    }

    /// Calculate probability of a certain allele.
    fn prob_allele<E>(
        &mut self,
        mut allele_params: E,
        edit_dist: &edit_distance::EditDistanceCalculation,
    ) -> LogProb
    where
        E: stats::pairhmm::EmissionParameters + pairhmm::RefBaseEmission,
    {
        let hit = edit_dist.calc_best_hit(&allele_params);
        if hit.dist() == 0 {
            // METHOD: In case of a perfect match, we just take the base quality product.
            // All alternative paths in the HMM will anyway be much worse.
            allele_params.read_emission().certainty_est()
        } else {
            // METHOD: We shrink the area to run the HMM against to an environment around the best
            // edit distance hits.
            allele_params.shrink_to_hit(&hit);

            // METHOD: Further, we run the HMM on a band around the best edit distance.
            self.pairhmm.prob_related(
                &self.gap_params,
                &allele_params,
                Some(hit.dist_upper_bound()),
            )
        }
    }
}

pub(crate) trait AltAlleleEmissionBuilder {
    type EmissionParams: stats::pairhmm::EmissionParameters + pairhmm::RefBaseEmission;

    fn build<'a>(
        &self,
        read_emission_params: &'a ReadEmission,
        ref_seq: &'a [u8],
    ) -> Self::EmissionParams;
}

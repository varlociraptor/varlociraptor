// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use anyhow::Result;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::bases::prob_read_base_miscall;
use crate::variants::evidence::observations::read_observation::Strand;
use crate::variants::evidence::realignment::Realignable;
use crate::variants::model;
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, Overlap, SingleEndEvidence, SingleLocus, Variant,
};

use super::ToVariantRepresentation;

#[derive(Debug)]
pub(crate) struct ImpreciseBreakendPair {
    breakends: Vec<ImpreciseBreakend>,
}

impl None {
    pub(crate) fn new(locus: genome::Locus, ref_base: u8) -> Self {
        None {
            locus: SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + 1,
            )),
            ref_base: ref_base.to_ascii_uppercase(),
        }
    }
}

impl Variant for None {
    type Evidence = SingleEndEvidence;
    type Loci = SingleLocus;

    fn is_valid_evidence(
        &self,
        evidence: &SingleEndEvidence,
        _: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        if let Overlap::Enclosing = self.locus.overlap(evidence, false) {
            Some(vec![0])
        } else {
            None
        }
    }

    fn loci(&self) -> &SingleLocus {
        &self.locus
    }

    fn allele_support(
        &self,
        read: &SingleEndEvidence,
        _: &AlignmentProperties,
        _: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        if let Some(qpos) = read
            .cigar_cached()
            .unwrap()
            // TODO expect u64 in read_pos
            .read_pos(self.locus.range().start as u32, false, false)?
        {
            let read_base = read.seq()[qpos as usize].to_ascii_uppercase();
            let base_qual = read.qual()[qpos as usize];
            let prob_miscall = prob_read_base_miscall(base_qual);

            let (prob_ref, prob_alt) = if read_base == self.ref_base {
                (prob_miscall.ln_one_minus_exp(), prob_miscall)
            } else {
                (prob_miscall, prob_miscall.ln_one_minus_exp())
            };

            let strand = if prob_ref != prob_alt {
                Strand::from_record_and_pos(read, qpos as usize)?
            } else {
                // METHOD: if record is not informative, we don't want to
                // retain its information (e.g. strand).
                Strand::None
            };

            Ok(Some(
                AlleleSupportBuilder::default()
                    .prob_ref_allele(prob_ref)
                    .prob_alt_allele(prob_alt)
                    .strand(strand)
                    .build()
                    .unwrap(),
            ))
        } else {
            // a read that spans a potential Ref site might have the respective position deleted (Cigar op 'D')
            // or reference skipped (Cigar op 'N'), and the library should not choke on those reads
            // but instead needs to know NOT to add those reads (as observations) further up
            Ok(None)
        }
    }

    fn prob_sample_alt(&self, _: &SingleEndEvidence, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

impl ToVariantRepresentation for None {
    fn to_variant_representation(&self) -> model::Variant {
        model::Variant::None
    }
}

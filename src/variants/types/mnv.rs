// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::cmp;
use std::ops::Range;
use std::rc::Rc;
use std::sync::Arc;

use anyhow::Result;

use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval};


use crate::default_ref_base_emission;
use crate::errors::Error;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils;
use crate::variants::evidence::bases::prob_read_base;
use crate::variants::evidence::observations::read_observation::Strand;
use crate::variants::evidence::realignment::edit_distance::is_explainable_by_error_rates;
use crate::variants::evidence::realignment::edit_distance::EditDistance;
use crate::variants::evidence::realignment::pairhmm::ReadEmission;
use crate::variants::evidence::realignment::pairhmm::RefBaseEmission;
use crate::variants::evidence::realignment::pairhmm::RefBaseVariantEmission;
use crate::variants::evidence::realignment::pairhmm::VariantEmission;
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::model;
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, Overlap, SingleEndEvidence, SingleLocus, Variant,
};

use super::ToVariantRepresentation;

#[derive(Debug)]
pub(crate) struct Mnv<R: Realigner> {
    locus: SingleLocus,
    ref_bases: Vec<u8>,
    alt_bases: Rc<Vec<u8>>,
    realigner: RefCell<R>,
}

impl<R: Realigner> Mnv<R> {
    pub(crate) fn new<L: genome::AbstractLocus>(
        locus: L,
        ref_bases: Vec<u8>,
        alt_bases: Vec<u8>,
        realigner: R,
    ) -> Self {
        Mnv {
            locus: SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + alt_bases.len() as u64,
            )),
            ref_bases: ref_bases.to_ascii_uppercase(),
            alt_bases: Rc::new(alt_bases.to_ascii_uppercase()),
            realigner: RefCell::new(realigner),
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.ref_bases.len()
    }
}

impl<R: Realigner> Realignable for Mnv<R> {
    fn alt_emission_params(
        &self,
        ref_buffer: Arc<reference::Buffer>,
        _: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<Box<dyn RefBaseVariantEmission>>> {
        let start = self.locus.range().start as usize;

        let ref_seq = ref_buffer.seq(self.locus.contig())?;

        let ref_seq_len = ref_seq.len();
        Ok(vec![Box::new(MnvEmissionParams {
            ref_seq,
            ref_offset: start.saturating_sub(ref_window),
            ref_end: cmp::min(start + self.alt_bases.len() + ref_window, ref_seq_len),
            alt_start: start,
            alt_end: self.locus.range().end as usize,
            alt_seq: Rc::clone(&self.alt_bases),
            ref_offset_override: None,
            ref_end_override: None,
        })])
    }
}

impl<R: Realigner> Variant for Mnv<R> {
    type Evidence = SingleEndEvidence;
    type Loci = SingleLocus;

    fn is_imprecise(&self) -> bool {
        false
    }

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
        alignment_properties: &AlignmentProperties,
        alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        if utils::contains_indel_op(&**read) || !alt_variants.is_empty() {
            // METHOD: reads containing indel operations should always be realigned,
            // as their support or non-support of the MNV might be an artifact
            // of the aligner. Also, if we have alt alignments here, we need to
            // realign as well since we need the multi-allelic case handling in the
            // realigner.
            Ok(Some(self.realigner.borrow_mut().allele_support(
                &**read,
                [&self.locus].iter(),
                self,
                alt_variants,
                alignment_properties,
            )?))
        } else {
            let mut prob_ref = LogProb::ln_one();
            let mut prob_alt = LogProb::ln_one();
            let mut prob_third = LogProb::ln_one();
            let aux_strand_info = utils::aux_tag_strand_info(read);
            let mut strand = Strand::None;
            let mut read_position = None;
            let mut alt_edit_dist = 0_u32;
            let read_emission = ReadEmission::new(read.seq(), read.qual(), None, None);
            let mut is_third_allele = false;

            for ((alt_base, ref_base), pos) in self
                .alt_bases
                .iter()
                .zip(self.ref_bases.iter())
                .zip(self.locus.range())
            {
                // TODO remove cast once read_pos uses u64
                if let Some(qpos) = read
                    .cigar_cached()
                    .unwrap()
                    .read_pos(pos as u32, false, false)?
                {
                    if read_position.is_none() {
                        // set first MNV position as read position
                        read_position = Some(qpos);
                    }
                    let read_base = unsafe { read.seq().decoded_base_unchecked(qpos as usize) };
                    let base_qual = unsafe { *read.qual().get_unchecked(qpos as usize) };

                    if read_base != *alt_base {
                        alt_edit_dist += 1;
                    }

                    let base_prob_alt = prob_read_base(read_base, *alt_base, base_qual);
                    let base_prob_ref = prob_read_base(read_base, *ref_base, base_qual);
                    let base_prob_third = prob_read_base(read_base, read_base, base_qual);

                    if base_prob_alt != base_prob_ref {
                        if let Some(strand_info) = aux_strand_info {
                            if let Some(s) = strand_info.get(qpos as usize) {
                                strand |= Strand::from_aux_item(*s)?;
                            } else {
                                return Err(Error::ReadPosOutOfBounds.into());
                            }
                        }
                    }

                    prob_ref += base_prob_ref;
                    prob_alt += base_prob_alt;
                    prob_third += base_prob_third;
                } else {
                    // a read that spans an SNV might have the respective position deleted (Cigar op 'D')
                    // or reference skipped (Cigar op 'N'), and the library should not choke on those reads
                    // but instead needs to know NOT to add those reads (as observations) further up
                    return Ok(None);
                }
            }

            if prob_alt > prob_ref {
                let explainable = is_explainable_by_error_rates(
                    alt_edit_dist as usize,
                    0,
                    0,
                    self.len(),
                    alignment_properties,
                    read_emission.error_rate(),
                );
                if alt_edit_dist > 0 && !explainable {
                    // METHOD: if the read supports ALT but has more substitutions compared to ALT
                    // than would be explainable by the expected error rate, it likely belongs to
                    // a third allele. In such a case, we override prob_ref by the probability that
                    // the read stems from an artificial third allele derived from its own sequence.
                    // This is a pragmatic, conservative approach to avoid false positives.
                    prob_ref = prob_third;
                    is_third_allele = true;
                }
            }

            if aux_strand_info.is_none() && prob_ref != prob_alt {
                // record global strand information
                // METHOD: if record is not informative, we don't want to
                // retain its information (e.g. strand).
                strand = Strand::from_record(read);
            }

            Ok(Some(
                AlleleSupportBuilder::default()
                    .strand(strand)
                    .prob_ref_allele(prob_ref)
                    .prob_alt_allele(prob_alt)
                    .read_position(read_position)
                    .third_allele_evidence(if is_third_allele {
                        Some(EditDistance(alt_edit_dist))
                    } else {
                        None
                    })
                    .build()
                    .unwrap(),
            ))
        }
    }

    fn prob_sample_alt(&self, _: &SingleEndEvidence, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

impl<R: Realigner> ToVariantRepresentation for Mnv<R> {
    fn to_variant_representation(&self) -> model::Variant {
        model::Variant::Mnv(self.alt_bases.to_vec())
    }
}

/// Emission parameters for PairHMM over insertion allele.
pub(crate) struct MnvEmissionParams {
    ref_seq: Arc<Vec<u8>>,
    ref_offset: usize,
    ref_end: usize,
    alt_start: usize,
    alt_end: usize, // exclusive end
    alt_seq: Rc<Vec<u8>>,
    ref_offset_override: Option<usize>,
    ref_end_override: Option<usize>,
}

impl RefBaseEmission for MnvEmissionParams {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;

        if i_ < self.alt_start || i_ >= self.alt_end {
            self.ref_seq[i_]
        } else {
            self.alt_seq[i_ - self.alt_start]
        }
    }

    fn variant_homopolymer_ref_range(&self) -> Option<Range<u64>> {
        None
    }

    fn variant_ref_range(&self) -> Option<Range<u64>> {
        Some((self.alt_start as u64)..(self.alt_end as u64))
    }

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }

    default_ref_base_emission!();
}

impl VariantEmission for MnvEmissionParams {
    fn is_homopolymer_indel(&self) -> bool {
        false
    }
}

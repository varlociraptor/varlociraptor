// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::cmp;
use std::rc::Rc;
use std::sync::Arc;

use anyhow::Result;
use bio::stats::pairhmm::EmissionParameters;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval};

use crate::errors::Error;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils;
use crate::variants::evidence::bases::prob_read_base;
use crate::variants::evidence::observation::Strand;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, Overlap, SingleEndEvidence, SingleLocus, Variant,
};
use crate::{default_emission, default_ref_base_emission};

pub(crate) struct MNV<R: Realigner> {
    locus: SingleLocus,
    ref_bases: Vec<u8>,
    alt_bases: Rc<Vec<u8>>,
    realigner: RefCell<R>,
}

impl<R: Realigner> MNV<R> {
    pub(crate) fn new<L: genome::AbstractLocus>(
        locus: L,
        ref_bases: Vec<u8>,
        alt_bases: Vec<u8>,
        realigner: R,
    ) -> Self {
        MNV {
            locus: SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + alt_bases.len() as u64,
            )),
            ref_bases: ref_bases.to_ascii_uppercase(),
            alt_bases: Rc::new(alt_bases.to_ascii_uppercase()),
            realigner: RefCell::new(realigner),
        }
    }
}

impl<'a, R: Realigner> Realignable<'a> for MNV<R> {
    type EmissionParams = MNVEmissionParams<'a>;

    fn alt_emission_params(
        &self,
        read_emission_params: Rc<ReadEmission<'a>>,
        ref_buffer: Arc<reference::Buffer>,
        _: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<MNVEmissionParams<'a>>> {
        let start = self.locus.range().start as usize;

        let ref_seq = ref_buffer.seq(self.locus.contig())?;

        let ref_seq_len = ref_seq.len();
        Ok(vec![MNVEmissionParams {
            ref_seq,
            ref_offset: start.saturating_sub(ref_window),
            ref_end: cmp::min(start + self.alt_bases.len() + ref_window, ref_seq_len),
            alt_start: start,
            alt_end: self.locus.range().end as usize,
            alt_seq: Rc::clone(&self.alt_bases),
            read_emission: read_emission_params,
        }])
    }
}

impl<R: Realigner> Variant for MNV<R> {
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
    ) -> Result<Option<AlleleSupport>> {
        if utils::contains_indel_op(&**read) {
            // METHOD: reads containing indel operations should always be realigned,
            // as their support or non-support of the MNV might be an artifact
            // of the aligner.
            Ok(Some(self.realigner.borrow_mut().allele_support(
                &**read,
                [&self.locus].iter(),
                self,
            )?))
        } else {
            let mut prob_ref = LogProb::ln_one();
            let mut prob_alt = LogProb::ln_one();
            let aux_strand_info = utils::aux_tag_strand_info(read);
            let mut strand = Strand::None;
            let mut read_position = None;

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

                    // METHOD: instead of considering the actual REF base, we assume that REF is whatever
                    // base the read has at this position (if not the ALT base). This way, we avoid biased
                    // allele frequencies at sites with multiple alternative alleles.
                    // Note that this is an approximation. The real solution would be to have multiple allele
                    // frequency variables in the likelihood function, but that would be computationally
                    // more demanding (leading to a combinatorial explosion).
                    // However, the approximation is pretty accurate, because it will only matter for true
                    // multiallelic cases. Sequencing errors won't have a severe effect on the allele frequencies
                    // because they are too rare.
                    let non_alt_base = if read_base != *alt_base {
                        read_base
                    } else {
                        *ref_base
                    };

                    let base_prob_alt = prob_read_base(read_base, *alt_base, base_qual);
                    let base_prob_ref = prob_read_base(read_base, non_alt_base, base_qual);

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
                } else {
                    // a read that spans an SNV might have the respective position deleted (Cigar op 'D')
                    // or reference skipped (Cigar op 'N'), and the library should not choke on those reads
                    // but instead needs to know NOT to add those reads (as observations) further up
                    return Ok(None);
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
                    .build()
                    .unwrap(),
            ))
        }
    }

    fn prob_sample_alt(&self, _: &SingleEndEvidence, _: &AlignmentProperties) -> LogProb {
        LogProb::ln_one()
    }
}

/// Emission parameters for PairHMM over insertion allele.
pub(crate) struct MNVEmissionParams<'a> {
    ref_seq: Arc<Vec<u8>>,
    ref_offset: usize,
    ref_end: usize,
    alt_start: usize,
    alt_end: usize, // exclusive end
    alt_seq: Rc<Vec<u8>>,
    read_emission: Rc<ReadEmission<'a>>,
}

impl<'a> RefBaseEmission for MNVEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;

        if i_ < self.alt_start || i_ >= self.alt_end {
            self.ref_seq[i_]
        } else {
            self.alt_seq[i_ - self.alt_start]
        }
    }

    default_ref_base_emission!();
}

impl<'a> EmissionParameters for MNVEmissionParams<'a> {
    default_emission!();

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }
}

impl<'a> bio::stats::pairhmm::Emission for MNVEmissionParams<'a> {
    #[inline]
    fn emission_x(&self, i: usize) -> u8 {
        self.ref_base(i)
    }

    #[inline]
    fn emission_y(&self, j: usize) -> u8 {
        unsafe {
            self.read_emission
                .read_seq
                .decoded_base_unchecked(self.read_emission.project_j(j))
        }
    }
}

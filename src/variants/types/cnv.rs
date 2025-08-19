use std::rc::Rc;

use anyhow::Result;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval};
use rust_htslib::bam::Record;
use statrs::distribution::{Discrete, Poisson};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::model;
use crate::variants::types::breakends::{
    Breakend, BreakendGroup, BreakendGroupBuilder, ExtensionModification, Join, Side,
};
use crate::variants::types::{AlleleSupport, DepthVariant, Evidence, MultiLocus, ReadVariant};

use super::ToVariantRepresentation;

#[derive(Debug)]
pub(crate) struct Cnv<R: Realigner> {
    breakends: BreakendGroup<R>,
    len: u64,
}

impl<R: Realigner> Cnv<R> {
    pub(crate) fn new(interval: genome::Interval, realigner: R, chrom_seq: &[u8]) -> Self {
        let mut breakend_group_builder = BreakendGroupBuilder::new();
        breakend_group_builder.realigner(realigner);

        let get_ref_allele = |pos: u64| &chrom_seq[pos as usize..(pos + 1) as usize];
        let get_locus = |pos| genome::Locus::new(interval.contig().to_owned(), pos);

        // Duplication:
        let ref_allele = get_ref_allele(interval.range().start);
        breakend_group_builder.push_breakend(Breakend::from_operations(
            get_locus(interval.range().start),
            ref_allele,
            ref_allele.to_owned(),
            Join::new(
                genome::Locus::new(interval.contig().to_owned(), interval.range().end - 1),
                Side::LeftOfPos,
                ExtensionModification::None,
            ),
            false,
            b"u",
            b"w",
        ));

        let ref_allele = get_ref_allele(interval.range().end - 1);
        breakend_group_builder.push_breakend(Breakend::from_operations(
            get_locus(interval.range().end - 1),
            ref_allele,
            ref_allele.to_owned(),
            Join::new(
                genome::Locus::new(interval.contig().to_owned(), interval.range().start),
                Side::RightOfPos,
                ExtensionModification::None,
            ),
            true,
            b"w",
            b"u",
        ));

        // Deletion:
        let ref_allele = get_ref_allele(interval.range().start - 1);
        breakend_group_builder.push_breakend(Breakend::from_operations(
            get_locus(interval.range().start - 1),
            ref_allele,
            ref_allele.to_owned(),
            Join::new(
                genome::Locus::new(interval.contig().to_owned(), interval.range().end),
                Side::RightOfPos,
                ExtensionModification::None,
            ),
            true,
            b"a",
            b"b",
        ));

        let ref_allele = get_ref_allele(interval.range().end);
        breakend_group_builder.push_breakend(Breakend::from_operations(
            get_locus(interval.range().end),
            ref_allele,
            ref_allele.to_owned(),
            Join::new(
                genome::Locus::new(interval.contig().to_owned(), interval.range().start - 1),
                Side::LeftOfPos,
                ExtensionModification::None,
            ),
            false,
            b"b",
            b"a",
        ));

        Cnv {
            breakends: breakend_group_builder.build().unwrap(),
            len: interval.range().end - interval.range().start,
        }
    }

    fn read_starts_in_cnv(&self, read_pos: i64) -> Option<Vec<usize>> {
        let cnv_start = self.breakends.loci()[0].range().start as i64;
        let cnv_end = self.breakends.loci()[1].range().end as i64;

        if read_pos < cnv_start || read_pos > cnv_end {
            None
        } else {
            Some(vec![0])
        }
    }
}

impl<R: Realigner> ReadVariant for Cnv<R> {
    fn is_imprecise(&self) -> bool {
        false
    }

    fn is_valid_evidence(
        &self,
        evidence: &Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        self.breakends
            .is_valid_evidence(evidence, alignment_properties)
    }

    fn loci(&self) -> &MultiLocus {
        self.breakends.loci()
    }

    // If I understand it correctly, I can completely overtake this implementation from inversions, since it just uses the breakends defined in the struct as genomic regions.
    // Afterwards it compute the probabilities of prob_ref/ prob_alt by computing the minimum edit distance between the read and the prob_ref/ prob_alt alleles.
    fn allele_support(
        &self,
        evidence: &Evidence,
        alignment_properties: &AlignmentProperties,
        alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        let support =
            self.breakends
                .allele_support(evidence, alignment_properties, alt_variants)?;

        Ok(support)
    }

    fn prob_sample_alt(
        &self,
        evidence: &Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        self.breakends
            .prob_sample_alt(evidence, alignment_properties)
    }
}

impl<R: Realigner> DepthVariant for Cnv<R> {
    fn is_imprecise(&self) -> bool {
        false
    }

    fn is_valid_evidence(
        &self,
        evidence: &Evidence,
        _: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        // TODO: Understand and consider REF, understand interval borders (Are they included in the breakend?)
        match evidence {
            Evidence::SingleEndSequencingRead(read) => self.read_starts_in_cnv(read.inner.core.pos),
            Evidence::PairedEndSequencingRead { left, right } => self
                .read_starts_in_cnv(left.inner.core.pos)
                .or(self.read_starts_in_cnv(right.inner.core.pos)),
        }
        // self.breakends
        //     .is_valid_evidence(evidence, alignment_properties)
    }

    fn loci(&self) -> &MultiLocus {
        self.breakends.loci()
    }

    fn allele_support(
        &self,
        evidence: &[Evidence],
        _: &AlignmentProperties,
        _: &[Box<dyn Realignable>],
    ) -> Result<Option<Vec<f64>>> {
        let cnv_start = self.breakends.loci()[0].range().start;
        let cnv_end = self.breakends.loci()[1].range().end - 1;
        let mut cnv_read_depth = vec![0.0; self.len as usize - 1];
        for ev in evidence {
            let starts = match ev {
                Evidence::SingleEndSequencingRead(read) => {
                    vec![read.inner.core.pos as u64]
                }
                Evidence::PairedEndSequencingRead { left, right } => {
                    vec![left.inner.core.pos as u64, right.inner.core.pos as u64]
                }
            };

            for start in starts {
                if start >= cnv_start && start < cnv_end {
                    let idx = (start - cnv_start) as usize;
                    if idx < cnv_read_depth.len() {
                        cnv_read_depth[idx] += 1.0;
                    }
                }
            }
        }
        Ok(Some(cnv_read_depth))
    }
}

impl<R: Realigner> ToVariantRepresentation for Cnv<R> {
    fn to_variant_representation(&self) -> model::Variant {
        model::Variant::Cnv(self.len)
    }
}

impl<R: Realigner> Realignable for Cnv<R> {
    fn alt_emission_params(
        &self,
        ref_buffer: std::sync::Arc<crate::reference::Buffer>,
        ref_interval: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<Box<dyn crate::variants::evidence::realignment::pairhmm::RefBaseVariantEmission>>>
    {
        self.breakends
            .alt_emission_params(ref_buffer, ref_interval, ref_window)
    }
}

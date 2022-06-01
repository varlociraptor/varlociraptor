use std::iter;
use std::ops::Deref;

use anyhow::Result;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::model;
use crate::variants::types::breakends::{
    Breakend, BreakendGroup, BreakendGroupBuilder, ExtensionModification, Join, Side,
};
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, Variant};

pub(crate) struct Inversion<R: Realigner>(BreakendGroup<R>);

impl<R: Realigner> Deref for Inversion<R> {
    type Target = BreakendGroup<R>;

    fn deref(&self) -> &BreakendGroup<R> {
        &self.0
    }
}

impl<R: Realigner> Inversion<R> {
    pub(crate) fn new(interval: genome::Interval, realigner: R, chrom_seq: &[u8]) -> Self {
        let mut breakend_group_builder = BreakendGroupBuilder::new();
        breakend_group_builder.realigner(realigner);

        let get_ref_allele = |pos: u64| &chrom_seq[pos as usize..(pos + 1) as usize];
        let get_locus = |pos| genome::Locus::new(interval.contig().to_owned(), pos);

        // Encode inversion via breakends, see VCF spec.
        let ref_allele = get_ref_allele(interval.range().start - 1);
        breakend_group_builder.push_breakend(Breakend::from_operations(
            get_locus(interval.range().start - 1),
            ref_allele,
            ref_allele.to_owned(),
            Join::new(
                genome::Locus::new(interval.contig().to_owned(), interval.range().end - 1),
                Side::LeftOfPos,
                ExtensionModification::ReverseComplement,
            ),
            true,
            b"w",
            b"u",
        ));

        let ref_allele = get_ref_allele(interval.range().start);
        breakend_group_builder.push_breakend(Breakend::from_operations(
            get_locus(interval.range().start),
            ref_allele,
            ref_allele.to_owned(),
            Join::new(
                genome::Locus::new(interval.contig().to_owned(), interval.range().end),
                Side::RightOfPos,
                ExtensionModification::ReverseComplement,
            ),
            false,
            b"v",
            b"x",
        ));

        let ref_allele = get_ref_allele(interval.range().end - 1);
        breakend_group_builder.push_breakend(Breakend::from_operations(
            get_locus(interval.range().end - 1),
            ref_allele,
            ref_allele.to_owned(),
            Join::new(
                genome::Locus::new(interval.contig().to_owned(), interval.range().start - 1),
                Side::LeftOfPos,
                ExtensionModification::ReverseComplement,
            ),
            true,
            b"u",
            b"w",
        ));

        let ref_allele = get_ref_allele(interval.range().end);
        breakend_group_builder.push_breakend(Breakend::from_operations(
            get_locus(interval.range().end),
            ref_allele,
            ref_allele.to_owned(),
            Join::new(
                genome::Locus::new(interval.contig().to_owned(), interval.range().start),
                Side::RightOfPos,
                ExtensionModification::ReverseComplement,
            ),
            false,
            b"x",
            b"v",
        ));

        Inversion(breakend_group_builder.build())
    }
}

impl<R: Realigner> Variant for Inversion<R> {
    type Evidence = PairedEndEvidence;
    type Loci = MultiLocus;

    fn is_valid_evidence(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        (**self).is_valid_evidence(evidence, alignment_properties)
    }

    fn loci(&self) -> &Self::Loci {
        (**self).loci()
    }

    fn allele_support(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
        alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        let support = (**self).allele_support(evidence, alignment_properties, alt_variants)?;

        Ok(support)
    }

    fn prob_sample_alt(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        (**self).prob_sample_alt(evidence, alignment_properties)
    }

    fn to_variant_representation<'a>(&'a self) -> Box<dyn Iterator<Item = model::Variant> + 'a> {
        Box::new(iter::once(model::Variant::Inversion(self.len)))
    }
}

impl<R: Realigner> Realignable for Inversion<R> {
    fn alt_emission_params(
        &self,
        ref_buffer: std::sync::Arc<crate::reference::Buffer>,
        ref_interval: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<Box<dyn crate::variants::evidence::realignment::pairhmm::RefBaseVariantEmission>>>
    {
        (**self).alt_emission_params(ref_buffer, ref_interval, ref_window)
    }
}

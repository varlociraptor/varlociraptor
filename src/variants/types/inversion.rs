use anyhow::Result;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::realignment::Realigner;
use crate::variants::types::breakends::{
    Breakend, BreakendGroup, BreakendGroupBuilder, ExtensionModification, Join, Side,
};
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, Variant};

#[derive(Derefable)]
pub(crate) struct Inversion(#[deref] BreakendGroup);

impl Inversion {
    pub(crate) fn new(interval: genome::Interval, realigner: Realigner, chrom_seq: &[u8]) -> Self {
        let reference_buffer = realigner.ref_buffer();
        let mut breakend_group_builder = BreakendGroupBuilder::default();
        breakend_group_builder.set_realigner(realigner);

        let get_ref_allele = |pos: u64| &chrom_seq[pos as usize..(pos + 1) as usize];
        let get_locus = |pos| genome::Locus::new(interval.contig().to_owned(), pos);

        // Encode inversion via breakends, see VCF spec.
        let ref_allele = get_ref_allele(interval.range().start - 1);
        breakend_group_builder.push_breakend(Breakend::from_operations(
            get_locus(interval.range().start - 1),
            ref_allele,
            ref_allele.to_owned(),
            Join::new(
                genome::Locus::new(
                    interval.contig().to_owned(),
                    interval.range().end - 1,
                ),
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
                genome::Locus::new(
                    interval.contig().to_owned(),
                    interval.range().start - 1,
                ),
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

        let breakend_group = breakend_group_builder.build(reference_buffer).unwrap();

        Inversion(breakend_group)
    }
}

impl Variant for Inversion {
    type Evidence = PairedEndEvidence;
    type Loci = MultiLocus;

    fn is_valid_evidence(&self, evidence: &Self::Evidence) -> Option<Vec<usize>> {
        (**self).is_valid_evidence(evidence)
    }

    fn loci(&self) -> &Self::Loci {
        (**self).loci()
    }

    fn allele_support(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> Result<Option<AlleleSupport>> {
        let support = (**self).allele_support(evidence, alignment_properties)?;

        Ok(support)
    }

    fn prob_sample_alt(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        (**self).prob_sample_alt(evidence, alignment_properties)
    }
}

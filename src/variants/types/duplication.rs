use anyhow::Result;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::realignment::Realigner;
use crate::variants::types::breakends::{
    Breakend, BreakendGroup, ExtensionModification, Operation, Side,
};
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, Variant};

#[derive(Derefable)]
pub(crate) struct Duplication(#[deref] BreakendGroup);

impl Duplication {
    pub(crate) fn new(interval: genome::Interval, realigner: Realigner, chrom_seq: &[u8]) -> Self {
        let mut breakend_group = BreakendGroup::new(realigner);

        let get_ref_allele = |pos: u64| &chrom_seq[pos as usize..pos as usize + 1];
        let get_locus = |pos| genome::Locus::new(interval.contig().to_owned(), pos);

        // Encode duplication via breakends, see VCF spec.

        let ref_allele = get_ref_allele(interval.range().start);
        breakend_group.push(Breakend::from_operations(
            get_locus(interval.range().start),
            ref_allele,
            [
                Operation::Join {
                    locus: genome::Locus::new(
                        interval.contig().to_owned(),
                        interval.range().end - 1,
                    ),
                    side: Side::LeftOfPos,
                    extension_modification: ExtensionModification::None,
                },
                Operation::Replacement(ref_allele.to_owned()),
            ],
            b"u",
            b"w",
        ));

        let ref_allele = get_ref_allele(interval.range().end - 1);
        breakend_group.push(Breakend::from_operations(
            get_locus(interval.range().end - 1),
            ref_allele,
            [
                Operation::Replacement(ref_allele.to_owned()),
                Operation::Join {
                    locus: genome::Locus::new(interval.contig().to_owned(), interval.range().start),
                    side: Side::RightOfPos,
                    extension_modification: ExtensionModification::None,
                },
            ],
            b"w",
            b"u",
        ));

        Duplication(breakend_group)
    }
}

impl Variant for Duplication {
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

use anyhow::Result;
use bio::{
    alignment::{Alignment, AlignmentOperation},
    pattern_matching::myers::Myers,
};
use bio_types::genome::{self, AbstractLocus};
use itertools::Itertools;

use crate::{reference, variants::model::Variant};

/// Return true if variant deletes bases from the reference (even if it is encoded as a replacement).
pub(crate) fn homopolymer_indel_len(
    variant: &Variant,
    locus: &genome::Locus,
    reference_buffer: &reference::Buffer,
) -> Result<Option<i8>> {
    let seq = reference_buffer.seq(locus.contig())?;
    let rpos = locus.pos() as usize;
    match variant {
        Variant::Deletion(l) => {
            if is_homopolymer_seq(&seq[rpos..rpos + *l as usize]) && *l < 256 {
                Ok(Some(-(*l as i8)))
            } else {
                Ok(None)
            }
        }
        Variant::Insertion(insseq) => {
            if is_homopolymer_seq(&insseq)
                && ((rpos < seq.len()
                    && extend_homopolymer_stretch(insseq[0], &mut seq[rpos..].iter()) > 0)
                    || (rpos > 0
                        && extend_homopolymer_stretch(insseq[0], &mut seq[..rpos].iter().rev())
                            > 0))
                && insseq.len() < 256
            {
                Ok(Some(insseq.len() as i8))
            } else {
                Ok(None)
            }
        }
        Variant::Replacement {
            ref ref_allele,
            ref alt_allele,
        } => {
            let (pattern, text) = if ref_allele.len() > alt_allele.len() {
                // deletion
                (alt_allele, ref_allele)
            } else {
                // insertion
                (ref_allele, alt_allele)
            };
            if pattern.len() <= 64 {
                let mut myers = Myers::<u64>::new(pattern);
                let mut aln = Alignment::default();
                let mut matches = myers.find_all(text, pattern.len() as u8);
                let mut best_dist = None;
                let mut best_aln = None;

                while matches.next_alignment(&mut aln) {
                    if best_dist.is_none() || best_dist.unwrap() > aln.score {
                        best_dist = Some(aln.score);
                        best_aln = Some(aln.clone());
                    }
                }
                if let Some(best_aln) = best_aln {
                    Ok(
                        HomopolymerIndelOperation::extract(text, pattern, &best_aln.operations)
                            .map(|op| op.len),
                    )
                } else {
                    // Pattern too long, unlikely as hell that this is a homopolymer artifact, hence just ignore.
                    Ok(None)
                }
            } else {
                // no similarity, cannot be a homopolymer error
                Ok(None)
            }
        }
        _ => Ok(None),
    }
}

#[derive(Clone, Debug, CopyGetters)]
#[getset(get_copy = "pub(crate)")]
pub(crate) struct HomopolymerIndelOperation {
    len: i8,
    text_pos: usize,
}

impl HomopolymerIndelOperation {
    /// Extract the homopolymer indel operation if there is exactly one in the given pattern compared to the text.
    pub(crate) fn extract(
        text: &[u8],
        pattern: &[u8],
        alignment: &[AlignmentOperation],
    ) -> Option<Self> {
        let mut rpos = 0;
        let mut qpos = 0;
        let mut homopolymer_indel_len = None;
        let mut text_pos = 0;
        for (op, stretch) in &alignment.iter().group_by(|op| *op) {
            let len = stretch.count();
            match op {
                AlignmentOperation::Match => {
                    qpos += len;
                    rpos += len;
                }
                AlignmentOperation::Subst => {
                    qpos += len;
                    rpos += len;
                }
                AlignmentOperation::Del => {
                    if len < 256 && is_homopolymer_seq(&text[rpos..rpos + len]) {
                        if homopolymer_indel_len.is_none() {
                            homopolymer_indel_len = Some(-(len as i8));
                            text_pos = rpos;
                        } else {
                            // METHOD: more complex indel situation, not considered for homopolymer error handling.
                            return None;
                        }
                    }
                    rpos += len;
                }
                AlignmentOperation::Ins => {
                    if len <= 256
                        && is_homopolymer_seq(&pattern[qpos..qpos + len])
                        && ((rpos < text.len()
                            && extend_homopolymer_stretch(pattern[qpos], &mut text[rpos..].iter())
                                > 0)
                            || (rpos > 0
                                && extend_homopolymer_stretch(
                                    pattern[qpos],
                                    &mut text[..rpos].iter().rev(),
                                ) > 0))
                    {
                        if homopolymer_indel_len.is_none() {
                            homopolymer_indel_len = Some(len as i8);
                            text_pos = rpos;
                        } else {
                            // METHOD: more complex indel situation, not considered for homopolymer error handling.
                            return None;
                        }
                    }
                    qpos += len;
                }
                AlignmentOperation::Xclip(l) => {
                    assert_eq!(len, 1);
                    rpos += l;
                }
                AlignmentOperation::Yclip(l) => {
                    assert_eq!(len, 1);
                    qpos += l;
                }
            }
        }

        homopolymer_indel_len.map(|len| HomopolymerIndelOperation { len, text_pos })
    }
}

pub(crate) fn is_homopolymer_seq(seq: &[u8]) -> bool {
    seq[1..].iter().all(|c| *c == seq[0])
}

pub(crate) fn extend_homopolymer_stretch(base: u8, seq: &mut dyn Iterator<Item = &u8>) -> usize {
    seq.take_while(|c| **c == base).count()
}

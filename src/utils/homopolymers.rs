use std::collections::HashMap;

use bio::{
    alignment::{pairwise::Aligner, pairwise::Scoring, Alignment, AlignmentOperation},
    pattern_matching::myers::Myers,
    stats::{LogProb, Prob},
};

use itertools::Itertools;

use crate::{estimation::alignment_properties::AlignmentProperties, variants::types::Variant};

#[derive(Clone, Debug, CopyGetters)]
#[getset(get_copy = "pub(crate)")]
pub(crate) struct HomopolymerIndelOperation {
    len: i8,
    text_pos: usize,
    base: u8,
}

impl HomopolymerIndelOperation {
    pub(crate) fn from_text_and_pattern_global(text: &[u8], pattern: &[u8]) -> Option<Self> {
        let (text, pattern, reverse_direction) = if text.len() < pattern.len() {
            (pattern, text, true)
        } else {
            (text, pattern, false)
        };

        if pattern.len() <= 256 {
            let mut aligner = Aligner::with_scoring(Scoring::from_scores(-1, -1, 1, -1));
            let mut best_aln = aligner.global(pattern, text);

            let mut ret =
                HomopolymerIndelOperation::from_alignment(&text, &pattern, &best_aln.operations);
            if reverse_direction {
                if let Some(op) = ret.as_mut() {
                    op.len *= -1;
                }
            }
            ret
        } else {
            None
        }
    }

    /// Extract the homopolymer indel operation if there is exactly one in the given pattern compared to the text.
    pub(crate) fn from_alignment(
        text: &[u8],
        pattern: &[u8],
        alignment: &[AlignmentOperation],
    ) -> Option<Self> {
        let mut rpos = 0;
        let mut qpos = 0;
        let mut homopolymer_indel_len = None;
        let mut homopolymer_base = None;
        let mut text_pos = 0;

        let is_extendable_stretch = |rpos, base| {
            (rpos < text.len() && extend_homopolymer_stretch(base, &mut text[rpos..].iter()) > 0)
                || (rpos > 0
                    && extend_homopolymer_stretch(base, &mut text[..rpos].iter().rev()) > 0)
        };
        dbg!(alignment);
        dbg!(text);
        dbg!(pattern);

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
                    if len < 256
                        && is_homopolymer_seq(&text[rpos..rpos + len])
                        && is_extendable_stretch(rpos, text[rpos])
                    {
                        if homopolymer_indel_len.is_none() {
                            homopolymer_indel_len = Some(-(len as i8));
                            homopolymer_base = Some(text[rpos]);
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
                        && is_extendable_stretch(rpos, pattern[qpos])
                    {
                        if homopolymer_indel_len.is_none() {
                            homopolymer_indel_len = Some(len as i8);
                            homopolymer_base = Some(pattern[qpos]);
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

        homopolymer_indel_len.map(|len| HomopolymerIndelOperation {
            len,
            text_pos,
            base: homopolymer_base.unwrap(),
        })
    }
}

pub(crate) fn is_homopolymer_seq(seq: &[u8]) -> bool {
    let base = seq[0].to_ascii_uppercase();
    seq[1..].iter().all(|c| c.to_ascii_uppercase() == base)
}

pub(crate) fn extend_homopolymer_stretch(base: u8, seq: &mut dyn Iterator<Item = &u8>) -> usize {
    let base = base.to_ascii_uppercase();
    seq.take_while(|c| c.to_ascii_uppercase() == base).count()
}

#[derive(Debug, Clone, CopyGetters)]
#[getset(get_copy = "pub(crate)")]
pub(crate) struct HomopolymerErrorModel {
    prob_homopolymer_insertion: LogProb,
    prob_homopolymer_deletion: LogProb,
    variant_homopolymer_indel_len: i8,
}

impl HomopolymerErrorModel {
    pub(crate) fn new<V>(variant: &V, alignment_properties: &AlignmentProperties) -> Option<Self>
    where
        V: Variant,
    {
        if let Some(variant_homopolymer_indel_len) = variant.homopolymer_indel_len() {
            let prob_homopolymer_error = |condition: &dyn Fn(i8) -> bool| {
                LogProb::ln_sum_exp(
                    &alignment_properties
                        .wildtype_homopolymer_error_model
                        .iter()
                        .filter_map(|(item_len, prob)| {
                            if condition(*item_len) {
                                Some(LogProb::from(Prob(*prob)))
                            } else {
                                None
                            }
                        })
                        .collect_vec(),
                )
            };

            let prob_homopolymer_insertion = prob_homopolymer_error(&|item_len| item_len > 0);
            let prob_homopolymer_deletion = prob_homopolymer_error(&|item_len| item_len < 0);
            let prob_total = prob_homopolymer_insertion.ln_add_exp(prob_homopolymer_deletion);

            let model = Some(HomopolymerErrorModel {
                prob_homopolymer_insertion: prob_homopolymer_insertion - prob_total,
                prob_homopolymer_deletion: prob_homopolymer_deletion - prob_total,
                variant_homopolymer_indel_len,
            });

            model
        } else {
            None
        }
    }
}

use std::collections::HashMap;

use bio::{
    alignment::{Alignment, AlignmentOperation},
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
}

impl HomopolymerIndelOperation {
    pub(crate) fn from_text_and_pattern(text: &[u8], pattern: &[u8]) -> Option<Self> {
        let (text, pattern, reverse_direction) = if text.len() < pattern.len() {
            (pattern, text, true)
        } else {
            (text, pattern, false)
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
                let mut ret = HomopolymerIndelOperation::from_alignment(
                    &text[best_aln.xstart..],
                    pattern,
                    &best_aln.operations,
                );
                if reverse_direction {
                    if let Some(op) = ret.as_mut() {
                        op.len *= -1;
                    }
                }
                ret
            } else {
                // Pattern too long, unlikely as hell that this is a homopolymer artifact, hence just ignore.
                None
            }
        } else {
            // no similarity, cannot be a homopolymer error
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

pub(crate) struct HomopolymerErrorModel {
    wildtype: HashMap<i8, LogProb>,
    artifact: HashMap<i8, LogProb>,
}

impl HomopolymerErrorModel {
    pub(crate) fn new<V>(variant: &V, alignment_properties: &AlignmentProperties) -> Option<Self>
    where
        V: Variant,
    {
        if let Some(variant_indel_len) = variant.homopolymer_indel_len() {
            let normalize = |mut adjusted: HashMap<_, _>| {
                let prob_total =
                    LogProb::ln_sum_exp(&adjusted.values().cloned().collect::<Vec<_>>());
                for prob in adjusted.values_mut() {
                    *prob -= prob_total;
                }
                adjusted
            };

            let artifact_model = alignment_properties
                .wildtype_homopolymer_error_model
                .iter()
                .map(|(item_indel_len, prob)| {
                    // METHOD: do not consider any indel lens that are not in the same direction as the variant indel (they will be attributed to the ref allele
                    // and should therefore not be expected to be associated with the alt allele).
                    let prob = if *item_indel_len == 0
                        || (variant_indel_len < 0 && *item_indel_len > 0)
                        || (variant_indel_len > 0 && *item_indel_len < 0)
                    {
                        0.0
                    } else {
                        *prob
                    };
                    (
                        item_indel_len - variant_indel_len,
                        LogProb::from(Prob(prob)),
                    )
                })
                .collect();

            let wildtype_model = alignment_properties
                .wildtype_homopolymer_error_model
                .iter()
                .map(|(item_indel_len, prob)| {
                    // METHOD: do not consider any indel lens that are not in the same direction as the variant indel (they will be attributed to the ref allele
                    // and should therefore not be expected to be associated with the alt allele).
                    let prob = if (variant_indel_len < 0 && *item_indel_len > 0)
                        || (variant_indel_len > 0 && *item_indel_len < 0)
                    {
                        0.0
                    } else {
                        *prob
                    };
                    (*item_indel_len, LogProb::from(Prob(prob)))
                })
                .collect();

            Some(HomopolymerErrorModel {
                artifact: normalize(artifact_model),
                wildtype: normalize(wildtype_model),
            })
        } else {
            None
        }
    }

    pub(crate) fn prob_wildtype_homopolymer_error(&self, indel_len: i8) -> LogProb {
        self.wildtype
            .get(&indel_len)
            .cloned()
            .unwrap_or(LogProb::ln_zero())
    }

    pub(crate) fn prob_artifact_homopolymer_error(&self, indel_len: i8) -> LogProb {
        self.artifact
            .get(&indel_len)
            .cloned()
            .unwrap_or(LogProb::ln_zero())
    }
}

use std::collections::HashMap;

use bio::{
    alignment::{pairwise::Aligner, pairwise::Scoring, AlignmentOperation},
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
            let mut aligner = Aligner::with_scoring(Scoring::from_scores(-2, -1, 1, -1));
            let best_aln = aligner.global(pattern, text);

            if !is_single_indel(&best_aln.operations) {
                // METHOD: if replacement contains either multiple indels or additional substitutions
                // it is a complex variant which we do not consider for homopolymer errors.
                // As we basically test whether the indel is in phase with the others in such
                // cases, we should be safe in any case, as homopolymer errors are random and should
                // not occur in phase with each other or substitutions.
                return None;
            }

            let mut ret =
                HomopolymerIndelOperation::from_alignment(text, &pattern, &best_aln.operations);
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
            let min_len = if rpos < (text.len() - 1) && text[rpos] == base {
                0
            } else {
                1
            };
            (rpos < (text.len() - 1)
                && extend_homopolymer_stretch(base, &mut text[rpos + 1..].iter()) > min_len)
                || (rpos > 0
                    && extend_homopolymer_stretch(base, &mut text[..rpos].iter().rev()) > min_len)
        };

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

fn is_single_indel(alignment: &[AlignmentOperation]) -> bool {
    let op_blocks = alignment
        .iter()
        .group_by(|op| *op)
        .into_iter()
        .filter(|(op, _stretch)| match op {
            AlignmentOperation::Del | AlignmentOperation::Ins | AlignmentOperation::Subst => true,
            AlignmentOperation::Match => false,
            _ => unreachable!("bug: unexpected alignment operation"),
        })
        .count();

    op_blocks == 1
}

pub(crate) fn is_homopolymer_seq(seq: &[u8]) -> bool {
    let base = seq[0].to_ascii_uppercase();
    seq[1..].iter().all(|c| c.to_ascii_uppercase() == base)
}

pub(crate) fn is_homopolymer_iter(mut seq: impl Iterator<Item = u8>) -> bool {
    if let Some(first) = seq.next() {
        let base = first.to_ascii_uppercase();
        seq.all(|c| c.to_ascii_uppercase() == base)
    } else {
        // if the seq is empty, then it's also a homopolymer
        true
    }
}

pub(crate) fn extend_homopolymer_stretch(base: u8, seq: &mut dyn Iterator<Item = &u8>) -> usize {
    let base = base.to_ascii_uppercase();
    seq.take_while(|c| c.to_ascii_uppercase() == base).count()
}

#[derive(Debug, Clone, CopyGetters)]
pub(crate) struct HomopolymerErrorModel {
    error_model: HashMap<i8, LogProb>,
    #[getset(get_copy = "pub(crate)")]
    variant_homopolymer_indel_len: i8,
}

impl HomopolymerErrorModel {
    pub(crate) fn new<V>(variant: &V, alignment_properties: &AlignmentProperties) -> Option<Self>
    where
        V: Variant,
    {
        if let Some(variant_homopolymer_indel_len) = variant.homopolymer_indel_len() {
            let is_valid = |item_len| item_len != 0 && item_len <= i8::MAX as i16 && item_len >= i8::MIN as i16;

            let prob_total = LogProb::ln_sum_exp(
                &alignment_properties
                    .wildtype_homopolymer_error_model
                    .iter()
                    .filter_map(|(item_len, prob)| {
                        if is_valid(*item_len) {
                            Some(LogProb::from(Prob(*prob)))
                        } else {
                            None
                        }
                    })
                    .collect_vec(),
            );

            let error_model = alignment_properties
            .wildtype_homopolymer_error_model
            .iter()
            .filter_map(|(item_len, prob)| {
                if is_valid(*item_len) {
                    Some((*item_len as i8, LogProb::from(Prob(*prob)) - prob_total))
                } else {
                    None
                }
            }).collect();

            Some(HomopolymerErrorModel {
                error_model,
                variant_homopolymer_indel_len,
            })
        } else {
            None
        }
    }

    pub(crate) fn prob_observable(&self, homopolymer_indel_len: i8) -> LogProb {
        assert!(homopolymer_indel_len != 0);
        self.error_model
            .get(&homopolymer_indel_len)
            .copied()
            .unwrap_or(LogProb::ln_zero())
    }
}

#[cfg(test)]
mod test_homopolymer {
    use crate::utils::homopolymers::HomopolymerIndelOperation;
    use bio_types::alignment::AlignmentOperation::*;
    #[test]
    fn test_homopolymer_indel_operation() {
        // Insertion after identical base
        let test_alignment = &[Match, Match, Ins, Match, Match];
        let test_homopolymer_indel =
            HomopolymerIndelOperation::from_alignment(b"ACGT", b"ACCGT", test_alignment);
        if test_homopolymer_indel.is_some() {
            panic!("Invalid homopolymer error detected");
        }

        // Insertion before identical base
        let test_alignment = &[Match, Ins, Match, Match, Match];
        let test_homopolymer_indel =
            HomopolymerIndelOperation::from_alignment(b"ACGT", b"ACCGT", test_alignment);
        if test_homopolymer_indel.is_some() {
            panic!("Invalid homopolymer error detected");
        }

        // Insertion at the beginning of the homopolymer
        let test_alignment = &[Match, Ins, Match, Match, Match];
        let test_homopolymer_indel =
            HomopolymerIndelOperation::from_alignment(b"GTTA", b"GTTTA", test_alignment);
        if test_homopolymer_indel.is_none() {
            panic!("Missed homopolymer error");
        }

        // Insertion in the middle of the homopolymer
        let test_alignment = &[Match, Match, Ins, Match, Match];
        let test_homopolymer_indel =
            HomopolymerIndelOperation::from_alignment(b"GTTA", b"GTTTA", test_alignment);
        if test_homopolymer_indel.is_none() {
            panic!("Missed homopolymer error ");
        }
        // Insertion at the end of the homopolymer
        let test_alignment = &[Match, Match, Match, Ins, Match];
        let test_homopolymer_indel =
            HomopolymerIndelOperation::from_alignment(b"GTTA", b"GTTTA", test_alignment);
        if test_homopolymer_indel.is_none() {
            panic!("Missed homopolymer error");
        }

        // Insertion of identical base at end of text
        let test_alignment = &[Match, Match, Match, Match, Ins];
        let test_homopolymer_indel =
            HomopolymerIndelOperation::from_alignment(b"ACGT", b"ACGTT", test_alignment);
        if test_homopolymer_indel.is_some() {
            panic!("Invalid homopolymer error detected");
        }
    }
}

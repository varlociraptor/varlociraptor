// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use bio::stats::{LogProb, PHREDProb, Prob};
use rust_htslib::bam::Record;
use std::rc::Rc;

lazy_static! {
    static ref PROB_CONFUSION: LogProb = LogProb::from(Prob(0.3333));
    pub(crate) static ref PROB_ANY: LogProb = LogProb::from(Prob(0.25));
}

/// Calculate probability of read_base given ref_base.
pub(crate) fn prob_read_base(mut read_base: u8, ref_base: u8, base_qual: u8) -> LogProb {
    read_base = read_base.to_ascii_uppercase();
    let ref_base = ref_base.to_ascii_uppercase();
    if read_base == ref_base {
        unsafe { *BASEQUAL_TO_PROB_CALL.get_unchecked(base_qual as usize) }
    // TODO: Since we use this method for the realignment HMM, we have special cases for N and IUPAC codes. This is new and we need to think about it if this causes issues.
    } else if read_base == b'N' {
        // METHOD: N means there can be anything, assuming a flat probability here.
        *PROB_ANY
    } else if iupac_contains(read_base, ref_base) {
        // read_base is an IUPAC ambiguity code (introduced by
        // base conversio), and ref_base is one of the bases it can represent.
        unsafe { *BASEQUAL_TO_PROB_CALL.get_unchecked(base_qual as usize) }
    } else {
        let prob_miscall = prob_read_base_miscall(base_qual);
        // TODO replace the second term with technology specific confusion matrix
        prob_miscall + *PROB_CONFUSION
    }
}

/// Encode a set of bases as an IUPAC ambiguity code
pub(crate) fn bases_to_iupac(bases: &[u8]) -> u8 {
    let mut bases: Vec<u8> = bases.iter().map(|b| b.to_ascii_uppercase()).collect();
    bases.sort_unstable();
    bases.dedup();
    match bases.as_slice() {
        [b'A'] => b'A',
        [b'C'] => b'C',
        [b'G'] => b'G',
        [b'T'] => b'T',
        [b'A', b'C'] => b'M',
        [b'A', b'G'] => b'R',
        [b'A', b'T'] => b'W',
        [b'C', b'G'] => b'S',
        [b'C', b'T'] => b'Y',
        [b'G', b'T'] => b'K',
        [b'C', b'G', b'T'] => b'B',
        [b'A', b'G', b'T'] => b'D',
        [b'A', b'C', b'T'] => b'H',
        [b'A', b'C', b'G'] => b'V',
        [b'A', b'C', b'G', b'T'] => b'N',
        _ => panic!("bases_to_iupac: invalid bases: {:?}", bases),
    }
}

/// True if `base` is `N` or an IUPAC ambiguity code, i.e. it does not denote a single concrete base and hence must not be counted as a mismatch.
#[inline]
pub(crate) fn is_ambiguous_base(base: u8) -> bool {
    matches!(
        base,
        b'N' | b'M' | b'R' | b'W' | b'S' | b'Y' | b'K' | b'B' | b'D' | b'H' | b'V'
    )
}

/// True if `ref_base` (a concrete A/C/G/T) is among the bases represented by
/// the IUPAC ambiguity code `read_base`.
#[inline]
pub(crate) fn iupac_contains(read_base: u8, ref_base: u8) -> bool {
    matches!(
        (read_base, ref_base),
        (b'M', b'A')
            | (b'M', b'C')
            | (b'R', b'A')
            | (b'R', b'G')
            | (b'W', b'A')
            | (b'W', b'T')
            | (b'S', b'C')
            | (b'S', b'G')
            | (b'Y', b'C')
            | (b'Y', b'T')
            | (b'K', b'G')
            | (b'K', b'T')
            | (b'B', b'C')
            | (b'B', b'G')
            | (b'B', b'T')
            | (b'D', b'A')
            | (b'D', b'G')
            | (b'D', b'T')
            | (b'H', b'A')
            | (b'H', b'C')
            | (b'H', b'T')
            | (b'V', b'A')
            | (b'V', b'C')
            | (b'V', b'G')
    )
}

/// Unpack miscall probability of read_base.
pub(crate) fn prob_read_base_miscall(base_qual: u8) -> LogProb {
    unsafe { *BASEQUAL_TO_PROB_MISCALL.get_unchecked(base_qual as usize) }
}

/// unpack miscall probability of read_base.
fn _prob_read_base_miscall(base_qual: u8) -> LogProb {
    LogProb::from(PHREDProb::from((base_qual) as f64))
}

lazy_static! {
    pub(crate) static ref BASEQUAL_TO_PROB_MISCALL: [LogProb; 256] = {
        let mut probs = [LogProb::ln_zero(); 256];
        for (qual, prob) in (0u8..=255u8).map(|qual| (qual, _prob_read_base_miscall(qual))) {
            probs[qual as usize] = prob;
        }
        probs
    };
    pub(crate) static ref BASEQUAL_TO_PROB_CALL: [LogProb; 256] = {
        let mut probs = [LogProb::ln_zero(); 256];
        for (qual, prob) in BASEQUAL_TO_PROB_MISCALL.iter().enumerate() {
            probs[qual] = prob.ln_one_minus_exp();
        }
        probs
    };
}

/// Returns the complement base for a given base
pub(crate) fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => base,
    }
}

/// Determines the orientation of a read based on its flags.
///
/// For single-end reads: returns true if the read is reverse-complemented.
/// For paired-end reads: returns true if the read is from the reverse strand
/// (either first-in-pair and reverse, or second-in-pair and forward).
///
/// # Arguments
/// * `read` - The sequencing read to check
///
/// # Returns
/// * `true` if the read is from the reverse strand, `false` otherwise
pub(crate) fn read_reverse_orientation(read: &Rc<Record>) -> bool {
    let read_paired = read.is_paired();
    let read_reverse = read.is_reverse();
    let read_first = read.is_first_in_template();
    if read_paired {
        read_reverse && read_first || !read_reverse && !read_first
    } else {
        read_reverse
    }
}

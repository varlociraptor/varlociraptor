// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use bio::stats::{LogProb, PHREDProb, Prob};

lazy_static! {
    static ref PROB_CONFUSION: LogProb = LogProb::from(Prob(0.3333));
}

/// Calculate probability of read_base given ref_base.
pub(crate) fn prob_read_base(read_base: u8, ref_base: u8, base_qual: u8) -> LogProb {
    if read_base.to_ascii_uppercase() == ref_base.to_ascii_uppercase() {
        unsafe { *BASEQUAL_TO_PROB_CALL.get_unchecked(base_qual as usize) }
    } else {
        let prob_miscall = prob_read_base_miscall(base_qual);
        // TODO replace the second term with technology specific confusion matrix
        prob_miscall + *PROB_CONFUSION
    }
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

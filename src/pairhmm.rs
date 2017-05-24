use std::mem;

use itertools::Itertools;

use bio::stats::LogProb;


/// A pair Hidden Markov Model for comparing sequences x and y as described by
/// Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). Biological Sequence Analysis.
/// Current Topics in Genome Analysis 2008. http://doi.org/10.1017/CBO9780511790492.
///
pub struct PairHMM {
    fm: [Vec<LogProb>; 2],
    fx: [Vec<LogProb>; 2],
    fy: [Vec<LogProb>; 2],
    prob_cols: Vec<LogProb>
}


impl PairHMM {
    pub fn new() -> Self {
        PairHMM {
            fm: [Vec::new(), Vec::new()],
            fx: [Vec::new(), Vec::new()],
            fy: [Vec::new(), Vec::new()],
            prob_cols: Vec::new()
        }
    }

    /// Calculate the probability of sequence x being related to y via any alignment.
    pub fn prob_related<M, X, Y>(
        &mut self,
        prob_gap_x: LogProb,
        prob_gap_y: LogProb,
        prob_gap_x_extend: LogProb,
        prob_gap_y_extend: LogProb,
        prob_emit_xy: &M,
        prob_emit_x: &X,
        prob_emit_y: &Y,
        len_x: usize,
        len_y: usize,
        free_start_gap_x: bool,
        free_end_gap_x: bool
    ) -> LogProb where
        M: Fn(usize, usize) -> LogProb,
        X: Fn(usize) -> LogProb,
        Y: Fn(usize) -> LogProb
    {
        println!("-----------");
        for k in 0..2 {
            self.fm[k].clear();
            self.fx[k].clear();
            self.fy[k].clear();
            self.prob_cols.clear();

            self.fm[k].resize(len_y + 1, LogProb::ln_zero());
            self.fx[k].resize(len_y + 1, LogProb::ln_zero());
            self.fy[k].resize(len_y + 1, LogProb::ln_zero());

            if free_end_gap_x {
                let c = (len_x * 3).saturating_sub(self.prob_cols.capacity());
                self.prob_cols.reserve_exact(c);
            }
        }

        let prob_no_gap = prob_gap_x.ln_add_exp(prob_gap_y);
        let prob_no_gap_x_extend = prob_gap_x_extend.ln_one_minus_exp();
        let prob_no_gap_y_extend = prob_gap_y_extend.ln_one_minus_exp();

        let mut prev = 0;
        let mut curr = 1;
        self.fm[prev][0] = LogProb::ln_one();

        // iterate over x
        for i in 0..len_x {
            if i > 0 && free_start_gap_x {
                // allow alignment to start from offset in x (if prob_offset_x is set accordingly)
                self.fm[prev][0] = LogProb::ln_one();
            }

            // iterate over y
            for j in 0..len_y {
                let j_ = j + 1;

                // match or mismatch
                self.fm[curr][j_] = prob_emit_xy(i, j) + LogProb::ln_sum_exp(&[
                    // coming from state M
                    prob_no_gap.ln_one_minus_exp() + self.fm[prev][j_ - 1],
                    // coming from state X
                    prob_no_gap_x_extend + self.fx[prev][j_ - 1],
                    // coming from state Y
                    prob_no_gap_y_extend + self.fy[prev][j_ - 1]
                ]);

                // gap in y
                self.fx[curr][j_] = prob_emit_x(i) + (
                    // open gap
                    prob_gap_y + self.fm[prev][j_]
                ).ln_add_exp(
                    // extend gap
                    prob_gap_y_extend + self.fx[prev][j_]
                );

                // gap in x
                self.fy[curr][j_] = prob_emit_y(j) + (
                    // open gap
                    prob_gap_x + self.fm[curr][j_ - 1]
                ).ln_add_exp(
                    // extend gap
                    prob_gap_x_extend + self.fy[curr][j_ - 1]
                );

            }

            if free_end_gap_x {
                // Cache column probabilities or simply record the last probability.
                // We can put all of them in one array since we simply have to sum in the end.
                // This is also good for numerical stability.
                self.prob_cols.push(self.fm[curr].last().unwrap().clone());
                self.prob_cols.push(self.fx[curr].last().unwrap().clone());
                self.prob_cols.push(self.fy[curr].last().unwrap().clone());
            }

            // TODO remove
            for p in &self.fm[curr] {
                print!("{:.2} ", p.exp())
            }
            println!("");

            // next column
            mem::swap(&mut curr, &mut prev);
            // reset next column to zeros
            for v in self.fm[curr].iter_mut() {
                *v = LogProb::ln_zero();
            }
        }

        if free_end_gap_x {
            LogProb::ln_sum_exp(&self.prob_cols)
        } else {
            LogProb::ln_sum_exp(&[
                self.fm[prev].last().unwrap().clone(),
                self.fx[prev].last().unwrap().clone(),
                self.fy[prev].last().unwrap().clone()
            ])
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use bio::stats::{Prob, LogProb};
    use constants;


    fn prob_emit_xy(x: &[u8], y: &[u8], i: usize, j: usize) -> LogProb {
        if x[i] == y[j] {
            LogProb::from(Prob(1.0) - constants::PROB_ILLUMINA_SUBST)
        } else {
            LogProb::from(constants::PROB_ILLUMINA_SUBST / Prob(3.0))
        }
    }

    fn prob_emit_x_or_y() -> LogProb {
        LogProb::from(Prob(1.0) - constants::PROB_ILLUMINA_SUBST)
    }


    #[test]
    fn test_same() {
        let x = b"AGCTCGATCGATCGATC";
        let y = b"AGCTCGATCGATCGATC";

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(
            LogProb::from(constants::PROB_ILLUMINA_INS),
            LogProb::from(constants::PROB_ILLUMINA_DEL),
            LogProb::ln_zero(),
            LogProb::ln_zero(),
            &|i, j| prob_emit_xy(x, y, i, j),
            &|_| prob_emit_x_or_y(),
            &|_| prob_emit_x_or_y(),
            x.len(),
            y.len(),
            false,
            false
        );
        assert!(*p <= 0.0);
        assert_relative_eq!(*p, 0.0, epsilon=0.1);
    }

    #[test]
    fn test_insertion() {
        let x = b"AGCTCGATCGATCGATC";
        let y = b"AGCTCGATCTGATCGATCT";

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(
            LogProb::from(constants::PROB_ILLUMINA_INS),
            LogProb::from(constants::PROB_ILLUMINA_DEL),
            LogProb::ln_zero(),
            LogProb::ln_zero(),
            &|i, j| prob_emit_xy(x, y, i, j),
            &|_| prob_emit_x_or_y(),
            &|_| prob_emit_x_or_y(),
            x.len(),
            y.len(),
            false,
            false
        );
        assert!(*p <= 0.0);
        assert_relative_eq!(p.exp(), constants::PROB_ILLUMINA_INS.powi(2), epsilon=1e-12);
    }

    #[test]
    fn test_deletion() {
        let x = b"AGCTCGATCTGATCGATCT";
        let y = b"AGCTCGATCGATCGATC";

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(
            LogProb::from(constants::PROB_ILLUMINA_INS),
            LogProb::from(constants::PROB_ILLUMINA_DEL),
            LogProb::ln_zero(),
            LogProb::ln_zero(),
            &|i, j| prob_emit_xy(x, y, i, j),
            &|_| prob_emit_x_or_y(),
            &|_| prob_emit_x_or_y(),
            x.len(),
            y.len(),
            false,
            false
        );
        assert!(*p <= 0.0);
        assert_relative_eq!(p.exp(), constants::PROB_ILLUMINA_DEL.powi(2), epsilon=1e-12);
    }

    #[test]
    fn test_mismatch() {
        let x = b"AGCTCGAGCGATCGATC";
        let y = b"TGCTCGATCGATCGATC";

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(
            LogProb::from(constants::PROB_ILLUMINA_INS),
            LogProb::from(constants::PROB_ILLUMINA_DEL),
            LogProb::ln_zero(),
            LogProb::ln_zero(),
            &|i, j| prob_emit_xy(x, y, i, j),
            &|_| prob_emit_x_or_y(),
            &|_| prob_emit_x_or_y(),
            x.len(),
            y.len(),
            false,
            false
        );
        assert!(*p <= 0.0);
        assert_relative_eq!(p.exp(), (constants::PROB_ILLUMINA_SUBST / Prob(3.0)).powi(2), epsilon=1e-7);
    }
}

use bio::stats::LogProb;


pub trait Parameters {
    /// Probability to start with an offset in x.
    fn prob_offset_x(&self) -> LogProb;

    /// Probability to open a gap in x.
    fn prob_gap_x(&self) -> LogProb;

    /// Probability to open a gap in y.
    fn prob_gap_y(&self) -> LogProb;

    /// Probability to extend a gap.
    fn prob_extend(&self) -> LogProb;

    /// Probability to emit x_i and y_j.
    fn prob_emit_xy(&self, i: usize, j: usize);

    /// Probability to emit x_i.
    fn prob_emit_x(&self, i: usize);

    /// Probability to emit y_j.
    fn prob_emit_y(&self, j: usize);

    fn len_x(&self) -> usize;

    fn len_y(&self) -> usize;
}



pub struct PairHMM<P: Parameters> {
    parameters: P,
    vm: [Vec<LogProb>; 2],
    vx: [Vec<LogProb>; 2],
    vy: [Vec<LogProb>; 2]
}


impl PairHMM {
    /// Calculate the probability of sequence y coming from sequence x.
    pub fn prob_any(&mut self) -> LogProb {
        for k in 0..2 {
            self.vm[k].clear();
            self.vx[k].clear();
            self.vy[k].clear();

            self.vm[k].resize(self.parameters.len_y() + 1, LogProb::ln_zero());
            self.vx[k].resize(self.parameters.len_y() + 1, LogProb::ln_zero());
            self.vy[k].resize(self.parameters.len_y() + 1, LogProb::ln_zero());
        }

        let prev = 0;
        let curr = 1;
        self.vm[prev][0] = LogProb::ln_one();
        for j in 1..self.parameters.len_x() {
            if j > 1 {
                // allow alignment to start from offset in x
                self.vm[prev][0] = self.parameters.prob_offset_x();
            }

            self.vm[curr][i]
        }
    }
}

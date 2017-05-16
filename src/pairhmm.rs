use bio::stats::LogProb;


pub trait Parameters {
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
}



pub struct PairHMM {
    
}

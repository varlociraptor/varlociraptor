use bio::stats::Prob;


// Single base insertion and deletion rates for R1 according to Schirmer et al.
// BMC Bioinformatics 2016, 10.1186/s12859-016-0976-y
pub static PROB_ILLUMINA_INS: Prob = Prob(2.8e-6);
pub static PROB_ILLUMINA_DEL: Prob = Prob(5.1e-6);
pub static PROB_ILLUMINA_SUBST: Prob = Prob(0.0021);

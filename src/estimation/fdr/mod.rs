use bio::stats::PHREDProb;

pub mod bh;
pub mod ev;

pub const ALPHAS: [f64; 32] = [
    0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04,
    0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8,
    0.9, 0.95, 1.0,
];

#[derive(Serialize)]
struct Record {
    alpha: f64,
    gamma: PHREDProb,
}

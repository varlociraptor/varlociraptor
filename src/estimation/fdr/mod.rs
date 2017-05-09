use bio::stats::PHREDProb;

pub mod bh;
pub mod ev;


pub const ALPHAS: [f64; 28] = [
    0.01,
    0.02,
    0.03,
    0.04,
    0.05,
    0.06,
    0.07,
    0.08,
    0.09,
    0.1,
    0.15,
    0.2,
    0.25,
    0.3,
    0.35,
    0.4,
    0.45,
    0.5,
    0.55,
    0.6,
    0.65,
    0.7,
    0.75,
    0.8,
    0.85,
    0.9,
    0.95,
    1.0
];


#[derive(RustcEncodable)]
struct Record {
    alpha: f64,
    gamma: PHREDProb
}

#[derive(Debug, new, Clone)]
pub(crate) struct DepthObservation {
    observed_depth: u64,
    expected_ref_depth: f64,
}

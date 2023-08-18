#[derive(new)]
pub(crate) struct Methylation {
    locus: SingleLocus,
}

impl Variant for Methylation {}

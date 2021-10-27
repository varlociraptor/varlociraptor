use std::collections::HashMap;

use bio::stats::LogProb;

#[derive(Debug)]
pub(crate) struct ArtifactParameters {
    homopolymer_error_model: HomopolymerErrorModel,
}

#[derive(Debug)]
pub(crate) struct HomopolymerErrorModel {
    none: HashMap<i8, LogProb>,
    some: HashMap<i8, LogProb>,
}

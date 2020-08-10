// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::types::Variant;

pub(crate) mod fragments;
pub(crate) mod reads;

pub(crate) use fragments::FragmentSamplingBias;
pub(crate) use reads::ReadSamplingBias;

pub(crate) trait SamplingBias: Variant {
    /// Number of bases that are feasible for overlapping the variant.
    fn feasible_bases(&self, read_len: u64, alignment_properties: &AlignmentProperties) -> u64;

    fn enclosable_len(&self) -> Option<u64>;
}

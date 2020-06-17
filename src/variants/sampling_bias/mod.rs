use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::Variant;

pub mod fragments;
pub mod reads;

pub use fragments::FragmentSamplingBias;
pub use reads::ReadSamplingBias;

pub trait SamplingBias: Variant {
    /// Number of bases that are feasible for overlapping the variant.
    fn feasible_bases(&self, read_len: u64, alignment_properties: &AlignmentProperties) -> u64 {
        if self.len() < alignment_properties.max_del_cigar_len as u64 {
            read_len
        } else {
            (read_len as f64 * alignment_properties.frac_max_softclip) as u64
        }
    }

    fn len(&self) -> u64;
}

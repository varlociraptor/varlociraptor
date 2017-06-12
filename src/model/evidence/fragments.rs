use std::cmp;

use model::Variant;


/// Calculate the number of positions a fragment can have in a given window according to
/// Sahlin et al. biorxiv 2015
/// (http://biorxiv.org/content/biorxiv/early/2015/08/04/023929.full.pdf).
///
/// # Arguments
///
/// * `insert_size` - observed insert size according to read mapper
/// * `left_read_len` - read length of left read
/// * `right_read_len` - read length of right read
/// * `window` - size of the considered sampling window (total, symmetric around the event)
pub fn n_fragment_positions(
    insert_size: u32,
    left_read_len: u32,
    right_read_len: u32,
    window: u32
) -> u32 {

    cmp::min(
        // The fragment needs to enclose the centerpoint on the reference.
        // TODO add +1 or not?
        cmp::max(insert_size as i32 - left_read_len as i32 - right_read_len as i32, 0),
        // The insert size needs to be shorter than the window on the reference sequence
        // (upper boundary)
        cmp::max(window as i32 - insert_size as i32 + 1, 0)
    ) as u32
}


#[cfg(test)]
mod tests {
    use super::*;

    fn _test_n_fragment_positions(insert_size: u32) -> u32 {
        // window to the left and right of the variant
        let window = 364;
        let read_len = 100;

        n_fragment_positions(insert_size, read_len, read_len, window)
    }

    #[test]
    fn test_n_fragment_positions_too_small() {
        let n = _test_n_fragment_positions(150);
        // Not enough space for placements around a centerpoint.
        assert_eq!(n, 0);

        let n = _test_n_fragment_positions(200);
        // Not enough space for placements around a centerpoint.
        assert_eq!(n, 0);
    }

    #[test]
    fn test_n_fragment_positions_exact() {
        let n = _test_n_fragment_positions(201);
        // Enough space for 1 placement around a centerpoint.
        assert_eq!(n, 1);
    }

    #[test]
    fn test_n_fragment_positions() {
        let n = _test_n_fragment_positions(202);
        // Enough space for 2 placements around a centerpoint.
        assert_eq!(n, 2);
    }

    #[test]
    fn test_n_fragment_positions_too_large() {
        let n = _test_n_fragment_positions(800);
        assert_eq!(n, 0);
    }
}


// use model::Variant;
// use pairhmm;
// use model::evidence;
//
//
// /// Calculate read evindence for an indel.
// pub struct IndelEvidence {
//     gap_params: evidence::reads::IndelGapParams,
//     pairhmm: pairhmm::PairHMM
// }
//
//
// impl IndelEvidence {
//     /// Create a new instance.
//     pub fn new(
//         prob_insertion_artifact: LogProb,
//         prob_deletion_artifact: LogProb,
//         prob_insertion_extend_artifact: LogProb,
//         prob_deletion_extend_artifact: LogProb,
//         window: u32
//     ) -> Self {
//         IndelEvidence {
//             gap_params: IndelGapParams {
//                 prob_insertion_artifact: prob_insertion_artifact,
//                 prob_deletion_artifact: prob_deletion_artifact,
//                 prob_insertion_extend_artifact: prob_insertion_extend_artifact,
//                 prob_deletion_extend_artifact: prob_deletion_extend_artifact
//             },
//             pairhmm: pairhmm::PairHMM::new(),
//             window: window
//         }
//     }
//
//     /// Calculate probability for reference and alternative allele.
//     pub fn prob(&mut self,
//         record: &bam::Record,
//         cigar: &CigarStringView,
//         start: u32,
//         variant: &Variant,
//         ref_seq: &[u8]
//     ) -> Result<(LogProb, LogProb), Box<Error>> {
//
//     }
// }

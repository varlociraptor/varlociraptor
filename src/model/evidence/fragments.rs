use std::cmp;
use std::error::Error;

use rgsl::randist::gaussian::ugaussian_P;
use bio::stats::LogProb;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::bam;

use model::Variant;
use model::sample::InsertSize;
use pairhmm;
use model::evidence;


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


/// Estimate insert size from read pair.
///
/// # Arguments
///
/// * `left` - left read of the pair
/// * `right` - right read of the pair
pub fn estimate_insert_size(left: &bam::Record, right: &bam::Record) -> u32 {
    let left_cigar = left.cigar();
    let right_cigar = right.cigar();
    let left_clips = evidence::Clips::trailing(&left_cigar);
    let right_clips = evidence::Clips::leading(&right_cigar);

    let left_end = (left_cigar.end_pos() as u32) + left_clips.both();
    let right_start = (right.pos() as u32).saturating_sub(
        right_clips.both()
    );
    // as defined by Torsten Seemann
    // (http://thegenomefactory.blogspot.nl/2013/08/paired-end-read-confusion-library.html)
    let inner_mate_distance = right_start as i32 - left_end as i32;
    println!("inner mate distance: {} {} {}", inner_mate_distance, right_start, left_end);

    let insert_size = inner_mate_distance +
                      evidence::read_len(left, &left_cigar) as i32 +
                      evidence::read_len(right, &right_cigar) as i32;
    assert!(insert_size > 0, "bug: insert size {} is smaller than zero", insert_size);

    insert_size as u32
}


/// Calculate read evindence for an indel.
pub struct IndelEvidence {
    gap_params: evidence::reads::IndelGapParams,
    pairhmm: pairhmm::PairHMM,
    window: u32,
    insert_size: InsertSize
}


impl IndelEvidence {
    /// Create a new instance.
    pub fn new(
        insert_size: InsertSize,
        prob_insertion_artifact: LogProb,
        prob_deletion_artifact: LogProb,
        prob_insertion_extend_artifact: LogProb,
        prob_deletion_extend_artifact: LogProb,
        window: u32
    ) -> Self {
        IndelEvidence {
            gap_params: evidence::reads::IndelGapParams {
                prob_insertion_artifact: prob_insertion_artifact,
                prob_deletion_artifact: prob_deletion_artifact,
                prob_insertion_extend_artifact: prob_insertion_extend_artifact,
                prob_deletion_extend_artifact: prob_deletion_extend_artifact
            },
            pairhmm: pairhmm::PairHMM::new(),
            window: window,
            insert_size: insert_size
        }
    }

    /// Calculate probability for reference and alternative allele.
    pub fn prob(&self,
        insert_size: u32,
        variant: &Variant
    ) -> Result<(LogProb, LogProb), Box<Error>> {
        let shift = match variant {
            &Variant::Deletion(_)  => variant.len() as f64,
            &Variant::Insertion(_) => -(variant.len() as f64),
            &Variant::SNV(_) => panic!("no fragment observations for SNV")
        };

        let p_alt = isize_pmf(
            insert_size as f64,
            self.insert_size.mean + shift,
            self.insert_size.sd
        );

        let p_ref = isize_pmf(
            insert_size as f64,
            self.insert_size.mean,
            self.insert_size.sd
        );

        Ok((p_ref, p_alt))
    }
}


/// as shown in http://www.milefoot.com/math/stat/pdfc-normaldisc.htm
pub fn isize_pmf(value: f64, mean: f64, sd: f64) -> LogProb {
    // TODO fix density in paper
    LogProb(
        (
            ugaussian_P((value + 0.5 - mean) / sd) -
            ugaussian_P((value - 0.5 - mean) / sd)
        ).ln()// - ugaussian_P(-mean / sd).ln()
    )
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

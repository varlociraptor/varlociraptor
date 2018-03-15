use std::cmp;
use std::error::Error;
use std::ops::Range;
use std::f64;

use itertools::Itertools;
use rgsl::randist::gaussian::ugaussian_P;
use bio::stats::LogProb;
use rust_htslib::bam;

use model::Variant;
use model::sample::InsertSize;
use model::evidence;
use model::evidence::observation::ProbSampleAlt;


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


/// Estimate the insert size from read pair projected on reference sequence including clips.
/// Note that this is is not the insert size of the real fragment but rather the insert size of
/// the alignment on the reference sequence.
///
/// # Arguments
///
/// * `left` - left read of the pair
/// * `right` - right read of the pair
pub fn estimate_insert_size(left: &bam::Record, right: &bam::Record) -> Result<u32, Box<Error>> {
    let left_cigar = left.cigar();
    let right_cigar = right.cigar();

    let aln = |rec: &bam::Record, cigar| -> Result<(u32, u32), Box<Error>> {
        Ok((
            (rec.pos() as u32).saturating_sub(evidence::Clips::leading(cigar).both()),
            cigar.end_pos()? as u32 + evidence::Clips::trailing(cigar).both()
        ))
    };

    let (left_start, left_end) = aln(left, &left_cigar)?;
    let (right_start, right_end) = aln(right, &right_cigar)?;
    // as defined by Torsten Seemann
    // (http://thegenomefactory.blogspot.nl/2013/08/paired-end-read-confusion-library.html)
    let inner_mate_distance = right_start as i32 - left_end as i32;
    debug!("inner mate distance: {} {} {}", inner_mate_distance, right_start, left_end);

    let insert_size = inner_mate_distance +
                      (left_end - left_start) as i32 +
                      (right_end - right_start) as i32;
    assert!(insert_size > 0, "bug: insert size {} is smaller than zero", insert_size);

    Ok(insert_size as u32)
}


/// Calculate read evindence for an indel.
pub struct IndelEvidence {
    insert_size: InsertSize
}


impl IndelEvidence {
    /// Create a new instance.
    pub fn new(
        insert_size: InsertSize
    ) -> Self {

        IndelEvidence {
            insert_size: insert_size
        }
    }

    /// Get range of insert sizes with probability above zero.
    /// We use 6 SDs around the mean.
    fn pmf_range(&self) -> Range<u32> {
        let m = self.insert_size.mean.round() as u32;
        let s = self.insert_size.sd.ceil() as u32 * 6;
        m.saturating_sub(s)..m + s
    }

    /// Get probability of given insert size from distribution shifted by the given value.
    fn pmf(&self,  insert_size: u32, shift: f64) -> LogProb {
        isize_pmf(
            insert_size as f64,
            self.insert_size.mean + shift,
            self.insert_size.sd
        )
    }

    /// Returns true if insert size is discriminative.
    pub fn is_discriminative(&self, variant: &Variant) -> bool {
        variant.len() as f64 > self.insert_size.sd
    }

    /// Calculate probability for reference and alternative allele.
    pub fn prob(&self,
        insert_size: u32,
        variant: &Variant
    ) -> Result<(LogProb, LogProb), Box<Error>> {
        let shift = match variant {
            &Variant::Deletion(_)  => variant.len() as f64,
            &Variant::Insertion(_) => {
                //(-(variant.len() as f64), variant.len())
                // We don't support insertions for now because it is not possible to reliably
                // detect that the fragment only overlaps the insertion at the inner read ends.
                // See Sample::overlap.
                panic!("bug: insert-size based probability for insertions is currently unsupported");
            },
            &Variant::SNV(_) => panic!("no fragment observations for SNV"),
            &Variant::None => panic!("no fragment observations for None")
        };

        let p_ref = self.pmf(insert_size, 0.0);
        let p_alt = self.pmf(insert_size, shift);

        Ok((p_ref, p_alt))
    }

    /// Probability to sample read from alt allele for each possible max softclip up to a given
    /// theoretical maximum.
    /// If variant is small enough to be in CIGAR, max_softclip should be set to None
    /// (i.e., ignored), and the method will only return one value.
    ///
    /// The key idea is to take the insert size distribution and calculate the expected probability
    /// by considering the number of valid placements over all placements for each possible insert
    /// size.
    pub fn prob_sample_alt(
        &self,
        left_read_len: u32,
        right_read_len: u32,
        variant: &Variant
    ) -> ProbSampleAlt {
        let expected_prob_enclose = |max_softclip, delta| {
            let read_offsets = left_read_len.saturating_sub(max_softclip) +
                               right_read_len.saturating_sub(max_softclip);

            let expected_p_alt = LogProb::ln_sum_exp(
                &self.pmf_range().filter_map(|x| {
                    if x < delta {
                        // if x is too small to enclose the variant, we skip it as it adds zero to the sum
                        None
                    } else {
                        Some(
                            self.pmf(x, 0.0) +
                            // probability to sample a valid placement
                            LogProb(
                                (x.saturating_sub(delta).saturating_sub(read_offsets) as f64).ln() -
                                (x.saturating_sub(delta) as f64).ln()
                            )
                        )
                    }
                }).collect_vec()
            );

            expected_p_alt
        };
        let expected_prob_overlap = |max_softclip, delta| {
            let n_alt_left = cmp::min(delta, left_read_len);
            let n_alt_right = cmp::min(delta, right_read_len);
            let n_alt_valid = cmp::min(n_alt_left, max_softclip) +
                              cmp::min(n_alt_right, max_softclip);

            LogProb((n_alt_valid as f64).ln() - ((n_alt_left + n_alt_right) as f64).ln())
        };

        let max_feasible = cmp::max(left_read_len, right_read_len);

        match variant {
            &Variant::Deletion(_)  => {
                // Deletion length does not affect sampling because the reads come from the allele
                // where the deleted sequence is not present ;-).
                let delta = 0;
                ProbSampleAlt::Dependent(
                    (0..max_feasible + 1).map(
                        |s| expected_prob_enclose(s, delta)
                    ).collect_vec()
                )
            },
            &Variant::Insertion(ref seq) => {
                let delta = seq.len() as u32;
                ProbSampleAlt::Dependent(
                    (0..max_feasible + 1).map(
                        |s| expected_prob_overlap(s, delta)
                    ).collect_vec()
                )
            },
            // for SNVs sampling is unbiased
            &Variant::SNV(_) | &Variant::None => return ProbSampleAlt::One
        }
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

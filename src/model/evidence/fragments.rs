use std::cmp;
use std::error::Error;
use std::f64;
use std::ops::Range;

use bio::stats::LogProb;
use itertools::Itertools;
use rgsl::randist::gaussian::ugaussian_P;
use rust_htslib::bam;

use estimation::alignment_properties::AlignmentProperties;
use model::evidence;
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
    window: u32,
) -> u32 {
    cmp::min(
        // The fragment needs to enclose the centerpoint on the reference.
        // TODO add +1 or not?
        cmp::max(
            insert_size as i32 - left_read_len as i32 - right_read_len as i32,
            0,
        ),
        // The insert size needs to be shorter than the window on the reference sequence
        // (upper boundary)
        cmp::max(window as i32 - insert_size as i32 + 1, 0),
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
    let left_cigar = left.cigar_cached().unwrap();
    let right_cigar = right.cigar_cached().unwrap();

    let aln = |rec: &bam::Record, cigar| -> Result<(u32, u32), Box<Error>> {
        Ok((
            (rec.pos() as u32).saturating_sub(evidence::Clips::leading(cigar).both()),
            cigar.end_pos()? as u32 + evidence::Clips::trailing(cigar).both(),
        ))
    };

    let (left_start, left_end) = aln(left, &left_cigar)?;
    let (right_start, right_end) = aln(right, &right_cigar)?;
    // as defined by Torsten Seemann
    // (http://thegenomefactory.blogspot.nl/2013/08/paired-end-read-confusion-library.html)
    let inner_mate_distance = right_start as i32 - left_end as i32;

    let insert_size =
        inner_mate_distance + (left_end - left_start) as i32 + (right_end - right_start) as i32;
    assert!(
        insert_size > 0,
        "bug: insert size {} is smaller than zero",
        insert_size
    );

    Ok(insert_size as u32)
}

/// Calculate read evindence for an indel.
pub struct IndelEvidence {}

impl IndelEvidence {
    /// Create a new instance.
    pub fn new() -> Self {
        IndelEvidence {}
    }

    /// Get range of insert sizes with probability above zero.
    /// We use 6 SDs around the mean.
    fn pmf_range(&self, alignment_properties: &AlignmentProperties) -> Range<u32> {
        let m = alignment_properties.insert_size().mean.round() as u32;
        let s = alignment_properties.insert_size().sd.ceil() as u32 * 6;
        m.saturating_sub(s)..m + s
    }

    /// Get probability of given insert size from distribution shifted by the given value.
    fn pmf(
        &self,
        insert_size: u32,
        shift: f64,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        isize_pmf(
            insert_size as f64,
            alignment_properties.insert_size().mean + shift,
            alignment_properties.insert_size().sd,
        )
    }

    /// Returns true if insert size is discriminative.
    pub fn is_discriminative(
        &self,
        variant: &Variant,
        alignment_properties: &AlignmentProperties,
    ) -> bool {
        variant.len() as f64 > alignment_properties.insert_size().sd
    }

    /// Calculate probability for reference and alternative allele.
    pub fn prob(
        &self,
        insert_size: u32,
        variant: &Variant,
        alignment_properties: &AlignmentProperties,
    ) -> Result<(LogProb, LogProb), Box<Error>> {
        let shift = match variant {
            &Variant::Deletion(_) => variant.len() as f64,
            &Variant::Insertion(_) => {
                //(-(variant.len() as f64), variant.len())
                // We don't support insertions for now because it is not possible to reliably
                // detect that the fragment only overlaps the insertion at the inner read ends.
                // See Sample::overlap.
                panic!(
                    "bug: insert-size based probability for insertions is currently unsupported"
                );
            }
            &Variant::SNV(_) => panic!("no fragment observations for SNV"),
            &Variant::None => panic!("no fragment observations for None"),
        };

        let p_ref = self.pmf(insert_size, 0.0, alignment_properties);
        let p_alt = self.pmf(insert_size, shift, alignment_properties);

        Ok((p_ref, p_alt))
    }

    /// Probability to sample read from alt allele for globally observed number of feasible
    /// positions.
    ///
    /// The key idea is to take the insert size distribution and calculate the expected probability
    /// by considering the number of valid placements over all placements for each possible insert
    /// size.
    pub fn prob_sample_alt(
        &self,
        left_read_len: u32,
        right_read_len: u32,
        variant: &Variant,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        // TODO for long reads always return one?
        let expected_prob = |left_feasible, right_feasible, delta_ref, delta_alt, enclose_only: bool| {
            let mut infeasible_read_pos_left = left_read_len.saturating_sub(left_feasible);
            let mut infeasible_read_pos_right = right_read_len.saturating_sub(right_feasible);
            if !enclose_only {
                // If we don't need enclosing, the fragment may also overlap on the far right or left.
                // Hence, the are less infeasible positions (by the maximal feasible overlap).
                infeasible_read_pos_left = infeasible_read_pos_left.saturating_sub(left_feasible);
                infeasible_read_pos_right = infeasible_read_pos_left.saturating_sub(right_feasible);
            }
            let infeasible_read_pos = infeasible_read_pos_left + infeasible_read_pos_right;


            // Calculate over each possible true insert size instead of the concrete insert size
            // of this fragment. This way we get a global value. Otherwise, we would obtain a
            // bias for each fragment that is too small to enclose the variant,
            // e.g., a ref fragment.
            let expected_p_alt = LogProb::ln_sum_exp(
                &self
                    .pmf_range(alignment_properties)
                    .filter_map(|x| {
                        let internal_segment = x.saturating_sub(left_read_len)
                                                .saturating_sub(right_read_len);
                        // Number of posititions in the internal segment where variant may not start,
                        // because it would lead to zero overlap.
                        let infeasible_internal_pos_alt = (internal_segment + 1).saturating_sub(delta_alt);
                        let infeasible_pos_alt = infeasible_read_pos + infeasible_internal_pos_alt;
                        let infeasible_pos_ref = (internal_segment + 1).saturating_sub(delta_ref);

                        let valid_pos_alt = x.saturating_sub(delta_alt)
                                             .saturating_sub(infeasible_pos_alt);
                        let valid_pos_ref = x.saturating_sub(infeasible_pos_ref);

                        if x <= delta_alt || x <= delta_alt + infeasible_pos_alt || x <= infeasible_pos_ref {
                            // if x is too small to enclose the variant, we skip it as it adds zero to the sum
                            None
                        } else {
                            let p = self.pmf(x, 0.0, alignment_properties) +
                            // probability to sample a valid placement
                            LogProb(
                                (
                                    valid_pos_alt as f64
                                ).ln() -
                                (valid_pos_ref as f64).ln()
                            );

                            assert!(
                                p.is_valid(),
                                "bug: invalid probability sampling probability {:?}", p
                            );
                            Some(p)
                        }
                    }).collect_vec(),
            ).cap_numerical_overshoot(0.0001);

            expected_p_alt
        };

        let left_feasible = alignment_properties.feasible_bases(left_read_len, variant);
        let right_feasible = alignment_properties.feasible_bases(right_read_len, variant);

        match variant {
            &Variant::Deletion(delta_ref) => {
                // Deletion length does not affect sampling from alt allele because the reads come
                // from the allele where the deleted sequence is not present.
                let delta_alt = 0;
                expected_prob(left_feasible, right_feasible, delta_ref, delta_alt, true)
            }
            &Variant::Insertion(ref seq) => {
                let delta_ref = 0;
                let delta_alt = seq.len() as u32;
                expected_prob(left_feasible, right_feasible, delta_ref, delta_alt, false)
            }
            // for SNVs sampling is unbiased
            &Variant::SNV(_) | &Variant::None => LogProb::ln_one(),
        }
    }
}

/// as shown in http://www.milefoot.com/math/stat/pdfc-normaldisc.htm
pub fn isize_pmf(value: f64, mean: f64, sd: f64) -> LogProb {
    // TODO fix density in paper
    LogProb(
        (ugaussian_P((value + 0.5 - mean) / sd) - ugaussian_P((value - 0.5 - mean) / sd)).ln(),
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

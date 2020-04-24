use rust_htslib::bam;
use anyhow::Result;

use crate::model::evidence;

/// Estimate the insert size from read pair projected on reference sequence including clips.
/// Note that this is is not the insert size of the real fragment but rather the insert size of
/// the alignment on the reference sequence.
///
/// # Arguments
///
/// * `left` - left read of the pair
/// * `right` - right read of the pair
pub fn estimate_insert_size(left: &bam::Record, right: &bam::Record) -> Result<u64> {
    let left_cigar = left.cigar_cached().unwrap();
    let right_cigar = right.cigar_cached().unwrap();

    let aln = |rec: &bam::Record, cigar| -> Result<(u64, u64)> {
        Ok((
            (rec.pos() as u64).saturating_sub(evidence::Clips::leading(cigar).both() as u64),
            cigar.end_pos() as u64 + evidence::Clips::trailing(cigar).both() as u64,
        ))
    };

    let (left_start, left_end) = aln(left, &left_cigar)?;
    let (right_start, right_end) = aln(right, &right_cigar)?;
    // as defined by Torsten Seemann
    // (http://thegenomefactory.blogspot.nl/2013/08/paired-end-read-confusion-library.html)
    let inner_mate_distance = right_start as i64 - left_end as i64;

    let insert_size =
        inner_mate_distance + (left_end - left_start) as i64 + (right_end - right_start) as i64;
    assert!(
        insert_size > 0,
        "bug: insert size {} is smaller than zero",
        insert_size
    );

    Ok(insert_size as u64)
}
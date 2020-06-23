// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::f64;
use std::u32;

use anyhow::Result;
use itertools::Itertools;
use ordered_float::NotNan;
use rust_htslib::bam::{self, record::Cigar};
use statrs::statistics::{OrderStatistics, Statistics};

#[derive(Clone, Debug, Copy, Deserialize, Serialize)]
pub(crate) struct AlignmentProperties {
    insert_size: Option<InsertSize>,
    pub(crate) max_del_cigar_len: u32,
    pub(crate) max_ins_cigar_len: u32,
    pub(crate) frac_max_softclip: f64,
    pub(crate) max_read_len: u32,
}

impl AlignmentProperties {
    /// Update maximum observed cigar operation lengths. Return whether any D, I, S, or H operation
    /// was found in the cigar string.
    pub(crate) fn update_max_cigar_ops_len(
        &mut self,
        record: &bam::Record,
        allow_hardclips: bool,
    ) -> (bool, bool) {
        let norm = |j| NotNan::new(j as f64 / record.seq().len() as f64).unwrap();

        let mut is_regular = true;
        let mut has_soft_clip = false;
        for c in record.cigar().iter() {
            match c {
                &Cigar::SoftClip(l) => {
                    let s = norm(l);
                    self.frac_max_softclip =
                        *cmp::max(s, NotNan::new(self.frac_max_softclip).unwrap());
                    is_regular = false;
                    has_soft_clip = true;
                }
                &Cigar::Del(l) => {
                    self.max_del_cigar_len = cmp::max(l, self.max_del_cigar_len);
                    is_regular = false;
                }
                &Cigar::Ins(l) => {
                    self.max_ins_cigar_len = cmp::max(l, self.max_ins_cigar_len);
                    is_regular = false;
                }
                &Cigar::HardClip(_) if !allow_hardclips => {
                    is_regular = false;
                }
                _ => continue,
            }
        }

        (is_regular, has_soft_clip)
    }

    /// Estimate `AlignmentProperties` from first 10000 fragments of bam file.
    /// Only reads that are mapped, not duplicates and where quality checks passed are taken.
    pub(crate) fn estimate<R: bam::Read>(bam: &mut R, allow_hardclips: bool) -> Result<Self> {
        let mut properties = AlignmentProperties {
            insert_size: None,
            max_del_cigar_len: 0,
            max_ins_cigar_len: 0,
            frac_max_softclip: 0.0,
            max_read_len: 0,
        };

        let mut record = bam::Record::new();
        let mut tlens = Vec::new();
        let mut max_read_len = 0;
        let mut max_mapq = 0;
        let mut i = 0;
        let mut skipped = 0;
        let mut n_soft_clip = 0;
        let mut n_not_useable = 0;
        while i <= 10000 {
            if skipped >= 100000 {
                warn!(
                    "\nWARNING: Stopping alignment property estimation after skipping 100.000\n\
                     records and inspecting {} records. You should have another look\n\
                     at your reads.\n",
                    i
                );

                break;
            }
            if !bam.read(&mut record)? {
                break;
            }

            // Records to skip without updating max_cigar_ops_len AND without incrementing the
            // counter (to keep looking for 10000 useful records for the estimation)
            if record.mapq() == 0
                || record.is_duplicate()
                || record.is_quality_check_failed()
                || record.is_unmapped()
            {
                skipped += 1;
                continue;
            }

            max_mapq = cmp::max(max_mapq, record.mapq());
            max_read_len = cmp::max(max_read_len, record.seq().len());

            let (is_regular, has_soft_clip) =
                properties.update_max_cigar_ops_len(&record, allow_hardclips);

            // Records to skip after updating max_cigar_ops_len, BUT without incrementing the
            // counter (to keep looking for 10000 useful records for the estimation)
            if !record.is_paired()
                || !record.is_first_in_template()
                || !(record.tid() == record.mtid())
                || record.is_mate_unmapped()
            {
                skipped += 1;
                continue;
            }

            if !is_regular {
                n_not_useable += 1;
                if has_soft_clip {
                    n_soft_clip += 1;
                }
            } else {
                // record insert size
                tlens.push(record.insert_size().abs() as f64);
            }

            i += 1;
        }

        if tlens.len() == 0 {
            warn!(
                "\nFound no records to use for estimating the insert size. Will assume\n\
                single end sequencing data and calculate deletion probabilities without\n\
                considering the insert size.\n\
                \n\
                If your data should be paired end, please consider manually providing\n\
                --alignment-properties, e.g. computed with `samtools stats`. Also,\n\
                the following counts of unusable records might indicate a source of\n\
                this problem:\n\n\
                - I, D, S or H CIGAR operation: {nu}\n\
                - S CIGAR (soft clip, e.g. due to UMIs or adapters): {sc}\n\
                \n\
                In addition, {nr} records were skipped in the estimation for one\n\
                of the following reasons:\n\
                - not paired\n\
                - not the first segment with regard to the template sequence\n\
                - mapping quality of 0\n\
                - marked as a duplicate\n\
                - mate mapped to different template (e.g. different chromosome)\n\
                - failed some quality check according to the 512 SAM flag\n\
                - mate unmapped\n\
                - record unmapped\n",
                nu = n_not_useable,
                sc = n_soft_clip,
                nr = skipped
            );
            properties.insert_size = None;
            Ok(properties)
        } else {
            let upper = tlens.percentile(95);
            let lower = tlens.percentile(5);
            let valid = tlens
                .into_iter()
                .filter(|l| *l <= upper && *l >= lower)
                .collect_vec();

            properties.insert_size = Some(InsertSize {
                mean: valid.iter().sum::<f64>() / valid.len() as f64,
                sd: valid.iter().std_dev(),
            });
            Ok(properties)
        }
    }

    pub(crate) fn insert_size(&self) -> &InsertSize {
        match &self.insert_size {
            Some(insert_size)  => insert_size,
            None => panic!("This is a bug.\n
                            AlignmentProperties has insert_size set to `None`, but you are trying to use it without checking for None.")
        }
    }
}

/// Expected insert size in terms of mean and standard deviation.
/// This should be estimated from unsorted(!) bam files to avoid positional biases.
#[derive(Copy, Clone, Debug, Default, Serialize, Deserialize)]
pub(crate) struct InsertSize {
    pub(crate) mean: f64,
    pub(crate) sd: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_estimate() {
        let mut bam = bam::Reader::from_path("tests/resources/tumor-first30000.bam").unwrap();

        let props = AlignmentProperties::estimate(&mut bam, false).unwrap();
        println!("{:?}", props);

        assert_relative_eq!(props.insert_size().mean, 311.9736111111111);
        assert_relative_eq!(props.insert_size().sd, 11.9001225301502);
        assert_eq!(props.max_del_cigar_len, 30);
        assert_eq!(props.max_ins_cigar_len, 12);
        assert_eq!(props.frac_max_softclip, 0.63);
    }

    #[test]
    fn test_estimate_all_reads_have_short_clips() {
        let mut bam =
            bam::Reader::from_path("tests/resources/tumor-first30000.reads_with_soft_clips.bam")
                .unwrap();

        let props = AlignmentProperties::estimate(&mut bam, false).unwrap();
        println!("{:?}", props);

        assert!(props.insert_size.is_none());
        assert_eq!(props.max_del_cigar_len, 2);
        assert_eq!(props.max_ins_cigar_len, 4);
        assert_eq!(props.frac_max_softclip, 0.63);
    }

    #[test]
    fn test_estimate_all_reads_single_end() {
        // this file contains only single-ended reads (artificially made single-ended with awk)
        let mut bam = bam::Reader::from_path(
            "tests/resources/tumor-first30000.bunch_of_reads_made_single_ended.bam",
        )
        .unwrap();

        let props = AlignmentProperties::estimate(&mut bam, false).unwrap();
        println!("{:?}", props);

        assert!(props.insert_size.is_none());
        assert_eq!(props.max_del_cigar_len, 0);
        assert_eq!(props.max_ins_cigar_len, 0);
        assert_eq!(props.frac_max_softclip, 0.03);
    }
}

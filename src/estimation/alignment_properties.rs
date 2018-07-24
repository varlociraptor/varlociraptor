use std::io;
use std::f64;
use std::u32;
use std::error::Error;
use std::str::FromStr;
use std::cmp;

use ordered_float::NotNaN;
use rust_htslib::bam::{self, record::Cigar};
use csv;
use statrs::statistics::{Statistics, OrderStatistics};
use itertools::Itertools;
use bio::stats::PHREDProb;

use model::Variant;


#[derive(Clone, Debug, Copy)]
pub struct AlignmentProperties {
    insert_size: InsertSize,
    max_del_cigar_len: u32,
    max_ins_cigar_len: u32,
    frac_max_softclip: f64,
    max_mapq: PHREDProb
}


impl AlignmentProperties {
    /// Constructs a dummy instance where all bases are feasible.
    pub fn default(insert_size: InsertSize) -> Self {
        AlignmentProperties {
            insert_size,
            max_del_cigar_len: 30,
            max_ins_cigar_len: 30,
            frac_max_softclip: 1.0,
            max_mapq: PHREDProb(60.0)
        }
    }

    /// Estimate `AlignmentProperties` from first 10000 fragments of bam file.
    /// Only reads that are mapped, not duplicates and where quality checks passed are taken.
    pub fn estimate<R: bam::Read>(bam: &mut R) -> Result<Self, Box<Error>> {
        let mut record = bam::Record::new();
        let mut tlens = Vec::new();
        let mut max_softclip = NotNaN::new(0.0).unwrap();
        let mut max_del_cigar_len = 0;
        let mut max_ins_cigar_len = 0;
        let mut max_mapq = 0;
        let mut i = 0;
        while i <= 30000 {
            if let Err(e) = bam.read(&mut record) {
                if e.is_eof() {
                    break;
                }
                return Err(Box::new(e));
            }
            if record.is_unmapped() || record.is_duplicate() || record.is_quality_check_failed() {
                continue
            }

            if !record.is_mate_unmapped() && record.is_first_in_template() && record.tid() == record.mtid() {
                tlens.push(record.insert_size().abs() as f64);
            }

            max_mapq = cmp::max(max_mapq, record.mapq());

            let norm = |j| NotNaN::new(j as f64 / record.seq().len() as f64).unwrap();

            for c in record.cigar().iter() {
                match c {
                    &Cigar::SoftClip(l) => {
                        let s = norm(l);
                        max_softclip = cmp::max(s, max_softclip);
                    },
                    &Cigar::Del(l) => {
                        max_del_cigar_len = cmp::max(l, max_del_cigar_len);
                    },
                    &Cigar::Ins(l) => {
                        max_ins_cigar_len = cmp::max(l, max_ins_cigar_len);
                    },
                    _ => continue
                }
            }
            i += 1;
        }

        let insert_size = {
            let upper = tlens.percentile(95);
            let lower = tlens.percentile(5);
            let mut valid = tlens.into_iter().filter(|l| *l <= upper && *l >= lower).collect_vec();
            InsertSize {
                mean: valid.median(),
                sd: valid.iter().std_dev()
            }
        };

        Ok(AlignmentProperties {
            insert_size,
            max_del_cigar_len,
            max_ins_cigar_len,
            frac_max_softclip: *max_softclip,
            max_mapq: PHREDProb(max_mapq as f64)
        })
    }

    // Number of bases that are feasible for overlapping the variant.
    pub fn feasible_bases(&self, read_len: u32, variant: &Variant) -> u32 {
        match variant {
            &Variant::Deletion(l) if l <= self.max_del_cigar_len => read_len,
            &Variant::Insertion(ref seq) if seq.len() as u32 <= self.max_ins_cigar_len => read_len,
            &Variant::SNV(_) => return read_len,
            &Variant::None => return read_len,
            _ => (read_len as f64 * self.frac_max_softclip) as u32
        }
    }

    pub fn insert_size(&self) -> &InsertSize {
        &self.insert_size
    }

    pub fn max_mapq(&self) -> PHREDProb {
        self.max_mapq
    }
}


/// Expected insert size in terms of mean and standard deviation.
/// This should be estimated from unsorted(!) bam files to avoid positional biases.
#[derive(Copy, Clone, Debug, Default)]
pub struct InsertSize {
    pub mean: f64,
    pub sd: f64
}


impl InsertSize {
    /// Obtain insert size from samtools stats output.
    pub fn from_samtools_stats<R: io::Read>(
        samtools_stats: &mut R
    ) -> Result<InsertSize, Box<Error>> {
        let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t')
                                               .comment(Some(b'#'))
                                               .has_headers(false)
                                               .flexible(true)
                                               .from_reader(samtools_stats);

        let mut insert_size = InsertSize::default();

        for rec in rdr.records() {
            let rec = rec?;
            if &rec[0] == "SN" {
                if &rec[1] == "insert size average:" {
                    insert_size.mean = f64::from_str(&rec[2])?;
                } else if &rec[1] == "insert size standard deviation:" {
                    insert_size.sd = f64::from_str(&rec[2])?;
                    break;
                }
            }
        }

        Ok(insert_size)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_estimate() {
        let mut bam = bam::Reader::from_path("tests/resources/tumor-first30000.bam").unwrap();

        let props = AlignmentProperties::estimate(&mut bam).unwrap();
        println!("{:?}", props);

        assert_relative_eq!(props.insert_size.mean, 312.0);
        assert_relative_eq!(props.insert_size.sd, 12.0);
    }


    #[test]
    fn test_parse_insert_size() {
        let insert_size = InsertSize::from_samtools_stats(
            &mut io::BufReader::new(
                fs::File::open("tests/resources/samtools_stats.example.txt").unwrap()
            )
        ).unwrap();
        assert_relative_eq!(insert_size.mean, 311.7);
        assert_relative_eq!(insert_size.sd, 15.5);
    }
}

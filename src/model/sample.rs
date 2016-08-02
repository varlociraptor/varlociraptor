use std::str;
use std::collections::HashMap;

use rgsl::randist::gaussian::gaussian_pdf;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Cigar;
use bio::stats::{logprobs, LogProb};


use model;
use model::Variant;


fn prob_mapping(mapq: u8) -> LogProb {
    logprobs::ln_1m_exp(logprobs::phred_to_log(mapq as f64))
}


/// An observation for or against a variant.
#[derive(Clone, Debug)]
pub struct Observation {
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    pub prob_mapping: LogProb,
    /// Probability that the read/read-pair comes from the alternative allele.
    pub prob_alt: LogProb,
    /// Probability that the read/read-pair comes from the reference allele.
    pub prob_ref: LogProb,
    /// Probability of the read/read-pair given that it has been mismapped.
    pub prob_mismapped: LogProb
}


/// Expected insert size in terms of mean and standard deviation.
/// This should be estimated from unsorted(!) bam files to avoid positional biases.
#[derive(Copy, Clone, Debug)]
pub struct InsertSize {
    pub mean: f64,
    pub sd: f64
}


/// A sequenced sample, e.g., a tumor or a normal sample.
pub struct Sample<P: model::priors::Model> {
    reader: bam::IndexedReader,
    pileup_window: u32,
    insert_size: InsertSize,
    prior_model: P
}


impl<P: model::priors::Model> Sample<P> {
    /// Create a new `Sample`.
    ///
    /// # Arguments
    ///
    /// * `bam` - BAM file with the aligned and deduplicated sequence reads.
    /// * `pileup_window` - Window around the variant that shall be search for evidence (e.g. 5000).
    /// * `insert_size` - estimated insert size
    /// * `prior_model` - Prior assumptions about allele frequency spectrum of this sample.
    pub fn new(bam: bam::IndexedReader, pileup_window: u32, insert_size: InsertSize, prior_model: P) -> Self {
        Sample {
            reader: bam,
            pileup_window: pileup_window,
            insert_size: insert_size,
            prior_model: prior_model
        }
    }

    /// Calculate prior probability for given allele frequency.
    pub fn prior_model(&self) -> &P {
        &self.prior_model
    }

    /// Extract observations for the given variant.
    pub fn extract_observations(&mut self, chrom: &[u8], start: u32, variant: Variant) -> Result<Vec<Observation>, String> {
        if let Some(tid) = self.reader.header.tid(chrom) {
            let mut observations = Vec::new();
            let (end, varpos) = match variant {
                Variant::Deletion(length)  => (start + length, (start + length / 2) as i32), // TODO do we really need two centerpoints?
                Variant::Insertion(length) => (start + length, start as i32)
            };
            let mut pairs = HashMap::new();

            if self.reader.seek(tid, start - self.pileup_window, end + self.pileup_window).is_err() {
                return Err(format!("Unable to seek to variant at {}:{}", str::from_utf8(chrom).unwrap(), start));
            }
            let mut record = bam::Record::new();
            loop {
                match self.reader.read(&mut record) {
                    Err(bam::ReadError::NoMoreRecord) => break,
                    Err(_) => return Err("Error reading BAM record.".to_owned()),
                    _ => ()
                }
                if record.is_duplicate() || record.is_unmapped() {
                    continue;
                }
                let pos = record.pos();
                let cigar = record.cigar();
                let end_pos = record.end_pos(&cigar);
                if pos <= varpos || end_pos >= varpos {
                    // overlapping alignment
                    observations.push(self.read_observation(&record, &cigar, varpos, variant));
                } else if end_pos <= varpos {
                    // need to check mate
                    // since the bam file is sorted by position, we can't see the mate first
                    if record.mpos() >= varpos {
                        pairs.insert(record.qname().to_owned(), record.mapq());
                    }
                } else if let Some(mate_mapq) = pairs.get(record.qname()) {
                    // mate already visited, and this read maps right of varpos
                    observations.push(self.fragment_observation(&record, *mate_mapq, variant));
                }
            }

            Ok(observations)
        } else {
            Err(format!("Sequence {} cannot be found in BAM", str::from_utf8(chrom).unwrap()))
        }
    }

    fn read_observation(&self, record: &bam::Record, cigar: &[Cigar], varpos: i32, variant: Variant) -> Observation {
        let (length, is_del) = match variant {
            Variant::Deletion(length) => (length, true),
            Variant::Insertion(length) => (length, false)
        };
        let is_close = |qpos: i32, qlength: u32| (varpos - (qpos + qlength as i32 / 2)).abs() <= 50;
        let is_similar_length = |qlength: u32| (length as i32 - qlength as i32).abs() <= 20;

        let mut qpos = record.pos();
        let prob_mapping = prob_mapping(record.mapq());
        let mut obs = Observation {
            prob_mapping: prob_mapping,
            prob_alt: 1.0f64.ln(),
            prob_ref: 0.0f64.ln(),
            prob_mismapped: 1.0 // if the read is mismapped, we assume sampling probability 1.0
        };
        for c in cigar {
            match c {
                &Cigar::Match(l) | &Cigar::RefSkip(l) | &Cigar::Equal(l) | &Cigar::Diff(l) => qpos += l as i32,
                &Cigar::Back(l) => qpos -= l as i32,
                &Cigar::Del(l) if is_del && is_similar_length(l) && is_close(qpos, l) => {
                    // supports alt allele
                    return obs;
                },
                &Cigar::Ins(l) if !is_del && is_similar_length(l) && is_close(qpos, l) => {
                    // supports alt allele
                    return obs;
                }
                _ => ()
            }
        }

        // support ref allele
        obs.prob_alt = 0.0f64.ln();
        obs.prob_ref = 1.0f64.ln();
        obs
    }

    fn fragment_observation(&self, record: &bam::Record, mate_mapq: u8, variant: Variant) -> Observation {
        let insert_size = record.insert_size();
        let shift = match variant {
            Variant::Deletion(length)  => length as f64,
            Variant::Insertion(length) => -(length as f64)
        };
        Observation {
            prob_mapping: prob_mapping(record.mapq()) + prob_mapping(mate_mapq),
            prob_alt: gaussian_pdf(insert_size as f64 - self.insert_size.mean, self.insert_size.sd),
            prob_ref: gaussian_pdf(insert_size as f64 - self.insert_size.mean + shift, self.insert_size.sd),
            prob_mismapped: 1.0 // if the fragment is mismapped, we assume sampling probability 1.0
        }
    }
}

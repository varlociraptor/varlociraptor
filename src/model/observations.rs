
use std::str;
use std::collections::HashMap;

use rgsl::randist::gaussian::gaussian_pdf;
use bio::stats::logprobs;
use bio::stats::LogProb;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Cigar;


fn prob_mapping(mapq: u8) -> LogProb {
    logprobs::ln_1m_exp(logprobs::phred_to_log(mapq as f64))
}


/// An observation for or against a variant.
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


/// Extract observations from BAM file.
pub struct BAMProcessor {
    reader: bam::IndexedReader,
    window: u32,
    exp_insert_size: f64,
    exp_insert_size_sd: f64
}


impl BAMProcessor {
    pub fn new(bam: bam::IndexedReader, window: u32, exp_insert_size: f64, exp_insert_size_sd: f64) -> Self {
        BAMProcessor {
            reader: bam,
            window: window,
            exp_insert_size: exp_insert_size,
            exp_insert_size_sd: exp_insert_size_sd
        }
    }

    pub fn extract_observations(&mut self, chrom: &[u8], start: u32, length: u32, is_del: bool) -> Result<Vec<Observation>, String> {
        if let Some(tid) = self.reader.header.tid(chrom) {
            let mut observations = Vec::new();
            let end = start + length;
            let varpos = if is_del { start + length / 2 } else { start } as i32; // TODO do we really need two centerpoints?
            let mut pairs = HashMap::new();

            if self.reader.seek(tid, start - self.window, end + self.window).is_err() {
                return Err(format!("Unable to seek to variant at {}:{}", str::from_utf8(chrom).unwrap(), start));
            }
            for record in self.reader.records() {
                if let Ok(record) = record {
                    let pos = record.pos();
                    let cigar = record.cigar();
                    let end_pos = record.end_pos(&cigar);
                    if pos < varpos && end_pos > varpos {
                        // overlapping alignment
                        observations.push(self.read_observation(&record, &cigar, start, length, is_del));
                    } else if end_pos < varpos {
                        // need to check mate
                        if record.mpos() >= varpos {
                            pairs.insert(record.qname().to_owned(), record.mapq());
                        }
                    } else if let Some(mate_mapq) = pairs.get(record.qname()) {
                        // mate already visited, and this read maps right of varpos
                        observations.push(self.fragment_observation(&record, *mate_mapq, length, is_del));
                    }
                } else {
                    return Err("Error reading BAM record.".to_owned());
                }
            }

            Ok(observations)
        } else {
            Err(format!("Sequence {} cannot be found in BAM", str::from_utf8(chrom).unwrap()))
        }
    }

    fn read_observation(&self, record: &bam::Record, cigar: &[Cigar], start: u32, length: u32, is_del: bool) -> Observation {
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
                &Cigar::Del(l) if is_del && qpos as u32 == start && l == length => {
                    // supports alt allele
                    return obs;
                },
                &Cigar::Ins(l) if !is_del && qpos as u32 == start && l == length => {
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

    fn fragment_observation(&self, record: &bam::Record, mate_mapq: u8, length: u32, is_del: bool) -> Observation {
        let insert_size = record.insert_size();
        let exp_var_insert_size = if is_del { self.exp_insert_size - length as f64 } else { self.exp_insert_size + length as f64 };
        Observation {
            prob_mapping: prob_mapping(record.mapq()) + prob_mapping(mate_mapq),
            prob_alt: gaussian_pdf(insert_size as f64 - self.exp_insert_size, self.exp_insert_size_sd),
            prob_ref: gaussian_pdf(insert_size as f64 - exp_var_insert_size, self.exp_insert_size_sd),
            prob_mismapped: 1.0 // if the fragment is mismapped, we assume sampling probability 1.0
        }
    }
}

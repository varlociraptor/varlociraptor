use std::str;
use std::collections::{HashMap, VecDeque, vec_deque};
use std::cmp;
use std::error::Error;
use std::f64::consts;

use rgsl::randist::gaussian::{gaussian_pdf, ugaussian_P};
use rgsl::error::erfc;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Cigar;
use bio::stats::{LogProb, PHREDProb, Prob};

use model;
use model::Variant;


pub fn prob_mapping(mapq: u8) -> LogProb {
    LogProb::from(PHREDProb(mapq as f64)).ln_one_minus_exp()
}


quick_error! {
    #[derive(Debug)]
    pub enum RecordBufferError {
        UnknownSequence(chrom: String) {
            description("unknown sequence")
            display("sequence {} cannot be found in BAM", chrom)
        }
    }
}


/// Ringbuffer of BAM records. This data structure ensures that no bam record is read twice while
/// extracting observations for given variants.
pub struct RecordBuffer {
    reader: bam::IndexedReader,
    inner: VecDeque<bam::Record>,
    window: u32,
    use_secondary: bool
}


unsafe impl Sync for RecordBuffer {}
unsafe impl Send for RecordBuffer {}


impl RecordBuffer {
    /// Create a new `RecordBuffer`.
    pub fn new(bam: bam::IndexedReader, window: u32, use_secondary: bool) -> Self {
        RecordBuffer {
            reader: bam,
            inner: VecDeque::with_capacity(window as usize * 2),
            window: window as u32,
            use_secondary: use_secondary
        }
    }

    /// Return end position of buffer.
    fn end(&self) -> Option<u32> {
        self.inner.back().map(|rec| rec.pos() as u32)
    }

    fn tid(&self) -> Option<i32> {
        self.inner.back().map(|rec| rec.tid())
    }

    /// Fill buffer around the given interval.
    pub fn fill(&mut self, chrom: &[u8], start: u32, end: u32) -> Result<(), Box<Error>> {
        if let Some(tid) = self.reader.header.tid(chrom) {
            let window_start = cmp::max(start as i32 - self.window as i32 - 1, 0) as u32;
            if self.inner.is_empty() || self.end().unwrap() < window_start || self.tid().unwrap() != tid as i32 {
                let end = self.reader.header.target_len(tid).unwrap();
                try!(self.reader.seek(tid, window_start, end));
                self.inner.clear();
            } else {
                // remove records too far left
                let to_remove = self.inner.iter().take_while(|rec| rec.pos() < window_start as i32).count();
                self.inner.drain(..to_remove);
            }

            // extend to the right
            for record in self.reader.records() {
                let record = try!(record);
                let pos = record.pos();
                if record.is_duplicate() || record.is_unmapped() {
                    continue;
                }
                if !self.use_secondary && record.is_secondary() {
                    continue;
                }
                self.inner.push_back(record);
                if pos > end as i32 + self.window as i32 {
                    break;
                }
            }

            Ok(())
        } else {
            Err(Box::new(RecordBufferError::UnknownSequence(str::from_utf8(chrom).unwrap().to_owned())))
        }
    }

    pub fn iter(&self) -> vec_deque::Iter<bam::Record> {
        self.inner.iter()
    }
}


/// An observation for or against a variant.
#[derive(Clone, Debug, RustcDecodable, RustcEncodable)]
pub struct Observation {
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    pub prob_mapping: LogProb,
    /// Probability that the read/read-pair comes from the alternative allele.
    pub prob_alt: LogProb,
    /// Probability that the read/read-pair comes from the reference allele.
    pub prob_ref: LogProb,
    /// Probability of the read/read-pair given that it has been mismapped.
    pub prob_mismapped: LogProb,
    /// Type of evidence.
    pub evidence: Evidence
}


#[derive(Clone, Debug, RustcDecodable, RustcEncodable)]
pub enum Evidence {
    /// Insert size of fragment
    InsertSize(u32),
    /// Alignment of a single read
    Alignment
}


/// Expected insert size in terms of mean and standard deviation.
/// This should be estimated from unsorted(!) bam files to avoid positional biases.
#[derive(Copy, Clone, Debug)]
pub struct InsertSize {
    pub mean: f64,
    pub sd: f64
}


/// A sequenced sample, e.g., a tumor or a normal sample.
pub struct Sample {
    record_buffer: RecordBuffer,
    use_fragment_evidence: bool,
    insert_size: InsertSize,
    likelihood_model: model::likelihood::LatentVariableModel,
    prob_spurious_isize: LogProb,
    prob_missed_insertion_alignment: LogProb,
    prob_missed_deletion_alignment: LogProb,
    prob_spurious_indel_alignment: LogProb,
    max_indel_dist: u32,
    max_indel_len_diff: u32
}


impl Sample {
    /// Create a new `Sample`.
    ///
    /// # Arguments
    ///
    /// * `bam` - BAM file with the aligned and deduplicated sequence reads.
    /// * `pileup_window` - Window around the variant that shall be search for evidence (e.g. 5000).
    /// * `use_fragment_evidence` - Whether to use read pairs that are left and right of variant.
    /// * `use_secondary` - Whether to use secondary alignments.
    /// * `insert_size` - estimated insert size
    /// * `prior_model` - Prior assumptions about allele frequency spectrum of this sample.
    /// * `likelihood_model` - Latent variable model to calculate likelihoods of given observations.
    /// * `prob_spurious_isize` - rate of wrongly reported insert size abberations (mapper dependent, BWA: 0.01332338, LASER: 0.05922201)
    /// * `prob_missed_insertion_alignment` - rate of missed insertion alignments if insertion is present (mapper dependent, BWA: 0.2138, LASER: 0.3460)
    /// * `prob_missed_deletion_alignment` - rate of missed deletion alignments if deletion is present (mapper dependent, BWA: 0.0310, LASER: 0.0964)
    /// * `prob_spurious_indel_alignment` - rate of wrongly reported indel alignments if indel is absent (mapper dependent)
    pub fn new(
        bam: bam::IndexedReader,
        pileup_window: u32,
        use_fragment_evidence: bool,
        use_secondary: bool,
        insert_size: InsertSize,
        likelihood_model: model::likelihood::LatentVariableModel,
        prob_spurious_isize: Prob,
        prob_missed_insertion_alignment: Prob,
        prob_missed_deletion_alignment: Prob,
        prob_spurious_indel_alignment: Prob
    ) -> Self {
        Sample {
            record_buffer: RecordBuffer::new(bam, pileup_window, use_secondary),
            use_fragment_evidence: use_fragment_evidence,
            insert_size: insert_size,
            likelihood_model: likelihood_model,
            prob_spurious_isize: LogProb::from(prob_spurious_isize),
            prob_missed_insertion_alignment: LogProb::from(prob_missed_insertion_alignment),
            prob_missed_deletion_alignment: LogProb::from(prob_missed_deletion_alignment),
            prob_spurious_indel_alignment: LogProb::from(prob_spurious_indel_alignment),
            max_indel_dist: 50,
            max_indel_len_diff: 20
        }
    }

    pub fn max_indel_dist(&mut self, dist: u32) -> &mut Self {
        self.max_indel_dist = dist;
        self
    }

    pub fn max_indel_len_diff(&mut self, diff: u32) -> &mut Self {
        self.max_indel_len_diff = diff;
        self
    }

    /// Return likelihood model.
    pub fn likelihood_model(&self) -> &model::likelihood::LatentVariableModel {
        &self.likelihood_model
    }

    /// Extract observations for the given variant.
    pub fn extract_observations(
        &mut self,
        chrom: &[u8],
        start: u32,
        variant: Variant
    ) -> Result<Vec<Observation>, Box<Error>> {
        let mut observations = Vec::new();
        let (end, varpos) = match variant {
            Variant::Deletion(length)  => (start + length, (start + length / 2) as i32), // TODO do we really need two centerpoints?
            Variant::Insertion(length) => (start + length, start as i32),
            Variant::SNV(_) => (start, start as i32)
        };
        let mut pairs = HashMap::new();
        let mut n_overlap = 0;

        // move window to the current variant
        debug!("Filling buffer...");
        try!(self.record_buffer.fill(chrom, start, end));
        debug!("Done.");

        // iterate over records
        for record in self.record_buffer.iter() {
            let pos = record.pos();
            let cigar = record.cigar();
            let end_pos = record.end_pos(&cigar);
            if pos < varpos && end_pos > varpos {
                // overlapping alignment
                observations.push(self.read_observation(&record, &cigar, varpos, variant));
                n_overlap += 1;
            } else if variant.has_fragment_evidence() &&
                      self.use_fragment_evidence &&
                      (record.is_first_in_template() || record.is_last_in_template())
            {
                if end_pos <= varpos {
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
        }
        debug!("Extracted observations ({} fragments, {} overlapping reads).", pairs.len(), n_overlap);
        debug!("{:?}", observations);
        Ok(observations)
    }

    fn read_observation(
        &self,
        record: &bam::Record,
        cigar: &[Cigar],
        varpos: i32,
        variant: Variant
    ) -> Observation {
        let is_close = |qpos: i32, qlength: u32| (varpos - (qpos + qlength as i32 / 2)).abs() as u32 <= self.max_indel_dist;
        let is_similar_length = |qlength: u32, varlength: u32| (varlength as i32 - qlength as i32).abs() as u32 <= self.max_indel_len_diff;
        let contains_varpos = |qpos: i32, qlength: u32|  qpos >= varpos && qpos + qlength as i32 <= varpos;
        let snv_prob_alt = |base: u8| {
            // position in read
            let i = (varpos - record.pos()) as usize;
            let prob_miscall = LogProb::from(PHREDProb::from(record.qual()[i] as f64));
            if record.seq()[i] == base {
                prob_miscall.ln_one_minus_exp()
            } else {
                prob_miscall
            }
        };

        let mut qpos = record.pos();
        let prob_mapping = prob_mapping(record.mapq());
        for c in cigar {
            match (c, variant) {
                // potential SNV evidence
                (&Cigar::Match(l), Variant::SNV(base)) |
                (&Cigar::Diff(l), Variant::SNV(base))
                if contains_varpos(qpos, l) => {
                    let p = snv_prob_alt(base);
                    return Observation {
                        prob_mapping: prob_mapping,
                        prob_alt: p,
                        prob_ref: p.ln_one_minus_exp(),
                        prob_mismapped: LogProb::ln_one(), // if the read is mismapped, we assume sampling probability 1.0
                        evidence: Evidence::Alignment
                    };
                },
                // potential indel evidence
                (&Cigar::Del(l), Variant::Deletion(length))
                if is_similar_length(l, length) && is_close(qpos, l) => {
                    return Observation {
                        prob_mapping: prob_mapping,
                        prob_alt: self.prob_missed_deletion_alignment.ln_one_minus_exp(),
                        prob_ref: self.prob_spurious_indel_alignment,
                        prob_mismapped: LogProb::ln_one(), // if the read is mismapped, we assume sampling probability 1.0
                        evidence: Evidence::Alignment
                    };
                },
                (&Cigar::Ins(l), Variant::Insertion(length))
                if is_similar_length(l, length) && is_close(qpos, l) => {
                    // supports alt allele
                    return Observation {
                        prob_mapping: prob_mapping,
                        prob_alt: self.prob_missed_insertion_alignment.ln_one_minus_exp(),
                        prob_ref: self.prob_spurious_indel_alignment,
                        prob_mismapped: LogProb::ln_one(), // if the read is mismapped, we assume sampling probability 1.0
                        evidence: Evidence::Alignment
                    };
                },
                // other
                (&Cigar::Match(l), _) |
                (&Cigar::RefSkip(l), _) |
                (&Cigar::Equal(l), _) |
                (&Cigar::Diff(l), _) |
                (&Cigar::Del(l), _) |
                (&Cigar::Ins(l), _) => qpos += l as i32,
                (&Cigar::Back(l), _) => qpos -= l as i32,
                _ => ()
            }
        }

        // support ref allele
        match variant {
            Variant::SNV(_) => {
                panic!("not yet implemented")
            },
            Variant::Deletion(_) => {
                Observation {
                    prob_mapping: prob_mapping,
                    prob_alt: self.prob_missed_deletion_alignment,
                    prob_ref: self.prob_spurious_indel_alignment.ln_one_minus_exp(),
                    prob_mismapped: LogProb::ln_one(), // if the read is mismapped, we assume sampling probability 1.0
                    evidence: Evidence::Alignment
                }
            },
            Variant::Insertion(_) => {
                Observation {
                    prob_mapping: prob_mapping,
                    prob_alt: self.prob_missed_insertion_alignment,
                    prob_ref: self.prob_spurious_indel_alignment.ln_one_minus_exp(),
                    prob_mismapped: LogProb::ln_one(), // if the read is mismapped, we assume sampling probability 1.0
                    evidence: Evidence::Alignment
                }
            }
        }
    }

    fn fragment_observation(
        &self,
        record: &bam::Record,
        mate_mapq: u8,
        variant: Variant
    ) -> Observation {
        let insert_size = record.insert_size().abs();
        let shift = match variant {
            Variant::Deletion(length)  => length as f64,
            Variant::Insertion(length) => -(length as f64),
            Variant::SNV(_) => panic!("no fragment observations for SNV")
        };
        let p_alt = (
            // case: correctly called indel
            self.prob_spurious_isize.ln_one_minus_exp() + isize_pmf(
                insert_size as f64,
                self.insert_size.mean + shift,
                self.insert_size.sd
            )
        ).ln_add_exp(
            // case: no indel, false positive call
            self.prob_spurious_isize +
            isize_pmf(
                insert_size as f64,
                self.insert_size.mean,
                self.insert_size.sd
            )
        );

        let obs = Observation {
            prob_mapping: prob_mapping(record.mapq()) + prob_mapping(mate_mapq),
            prob_alt: p_alt,
            prob_ref: isize_pmf(insert_size as f64, self.insert_size.mean, self.insert_size.sd),
            prob_mismapped: LogProb::ln_one(), // if the fragment is mismapped, we assume sampling probability 1.0
            evidence: Evidence::InsertSize(insert_size as u32)
        };

        obs
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


/// Continuous normal density (obsolete).
pub fn isize_density(value: f64, mean: f64, sd: f64) -> LogProb {
    LogProb(gaussian_pdf(value - mean, sd).ln())
}


/// Manual normal density (obsolete, we can use GSL (see above)).
pub fn isize_density_louis(value: f64, mean: f64, sd: f64) -> LogProb {
    let mut p = 0.5 / (1.0 - 0.5 * erfc((mean + 0.5) / sd * consts::FRAC_1_SQRT_2));
    p *= erfc((-value - 0.5 + mean)/sd * consts::FRAC_1_SQRT_2) - erfc((-value + 0.5 + mean)/sd * consts::FRAC_1_SQRT_2);

    LogProb(p.ln())
}


pub fn isize_mixture_density_louis(value: f64, d: f64, mean: f64, sd: f64, rate: f64) -> LogProb {
    let p = 0.5 / ( rate*(1.0 - 0.5*erfc((mean + 0.5)/sd*consts::FRAC_1_SQRT_2)) + (1.0 - rate)*(1.0 - 0.5*erfc((mean + d + 0.5)/sd* consts::FRAC_1_SQRT_2)) );
    LogProb((p * (
        rate*( erfc((-value - 0.5 + mean)/sd*consts::FRAC_1_SQRT_2) - erfc((-value + 0.5 + mean)/sd*consts::FRAC_1_SQRT_2) ) + (1.0 - rate)*( erfc((-value - 0.5 + mean + d)/sd*consts::FRAC_1_SQRT_2) - erfc((-value + 0.5 + mean + d)/sd*consts::FRAC_1_SQRT_2) )
    )).ln())
}

#[cfg(test)]
mod tests {
    use super::*;
    use model;
    use likelihood;

    use itertools::Itertools;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use bio::stats::{LogProb, PHREDProb, Prob};

    #[test]
    fn test_isize_density() {
        let d1 = isize_density(300.0, 312.0, 15.0);
        let d2 = isize_pmf(300.0, 312.0, 15.0);
        let d3 = isize_density_louis(300.0, 312.0, 15.0);
        println!("{} {} {}", *d1, *d2, *d3);

        let d_mix = isize_mixture_density_louis(212.0, -100.0, 312.0, 15.0, 0.05);
        let rate = LogProb(0.05f64.ln());
        let p_alt = (
            // case: correctly called indel
            rate.ln_one_minus_exp() + isize_pmf(
                212.0,
                312.0 - 100.0,
                15.0
            )
        ).ln_add_exp(
            // case: no indel, false positive call
            rate +
            isize_pmf(
                212.0,
                312.0,
                15.0
            )
        );

        println!("{} {}", d_mix.exp(), p_alt.exp());
    }

    #[test]
    fn test_prob_mapping() {
        assert_relative_eq!(prob_mapping(0).exp(), 0.0);
        assert_relative_eq!(prob_mapping(10).exp(), 0.9);
        assert_relative_eq!(prob_mapping(20).exp(), 0.99);
        assert_relative_eq!(prob_mapping(30).exp(), 0.999);
    }

    fn setup_sample(isize_mean: f64) -> Sample {
        Sample::new(
            bam::IndexedReader::new(&"tests/indels.bam").unwrap(),
            2500,
            true,
            true,
            InsertSize { mean: isize_mean, sd: 20.0 },
            likelihood::LatentVariableModel::new(1.0),
            Prob(0.0),
            Prob(0.0),
            Prob(0.0),
            Prob(0.0)
        )
    }

    #[test]
    fn test_read_observation_indel() {
        let variant = model::Variant::Insertion(10);
        // insertion starts at 528 + 20 and has length 10
        let varpos = 528 + 20 + 5;

        let sample = setup_sample(150.0);
        let bam = bam::Reader::new(&"tests/indels.bam").unwrap();
        let records = bam.records().collect_vec();

        let true_alt_probs = [LogProb::ln_one(), LogProb::ln_one(), LogProb::ln_zero(), LogProb::ln_one(), LogProb::ln_zero()];

        for (record, true_alt_prob) in records.into_iter().zip(true_alt_probs.into_iter()) {
            let record = record.unwrap();
            let cigar = record.cigar();
            let obs = sample.read_observation(&record, &cigar, varpos, variant);
            assert_relative_eq!(*obs.prob_alt, **true_alt_prob);
            assert_relative_eq!(*obs.prob_ref, *obs.prob_alt.ln_one_minus_exp());
            assert_relative_eq!(*obs.prob_mapping, *(LogProb::from(PHREDProb(60.0)).ln_one_minus_exp()));
            assert_relative_eq!(*obs.prob_mismapped, *LogProb::ln_one());
        }
    }

    #[test]
    fn test_fragment_observation_no_evidence() {
        let sample = setup_sample(150.0);
        let bam = bam::Reader::new(&"tests/indels.bam").unwrap();
        let records = bam.records().map(|rec| rec.unwrap()).collect_vec();

        for varlen in &[0, 5, 10, 100] {
            println!("varlen {}", varlen);
            println!("insertion");
            let variant = model::Variant::Insertion(*varlen);
            for record in &records {
                let obs = sample.fragment_observation(record, 60u8, variant);
                println!("{:?}", obs);
                if *varlen == 0 {
                    assert_relative_eq!(*obs.prob_ref, *obs.prob_alt);
                } else {
                    assert!(obs.prob_ref > obs.prob_alt);
                }
            }
            println!("deletion");
            let variant = model::Variant::Deletion(*varlen);
            for record in &records {
                let obs = sample.fragment_observation(record, 60u8, variant);
                println!("{:?}", obs);
                if *varlen == 0 {
                    assert_relative_eq!(*obs.prob_ref, *obs.prob_alt);
                } else {
                    assert!(obs.prob_ref > obs.prob_alt);
                }
            }
        }
    }

    #[test]
    fn test_fragment_observation_evidence() {
        let bam = bam::Reader::new(&"tests/indels.bam").unwrap();
        let records = bam.records().map(|rec| rec.unwrap()).collect_vec();

        println!("deletion");
        let sample = setup_sample(100.0);
        let variant = model::Variant::Deletion(50);
        for record in &records {
            let obs = sample.fragment_observation(record, 60u8, variant);
            println!("{:?}", obs);
            assert_relative_eq!(obs.prob_ref.exp(), 0.0, epsilon=0.001);
            assert!(obs.prob_alt > obs.prob_ref);
        }

        println!("insertion");
        let sample = setup_sample(200.0);
        let variant = model::Variant::Insertion(50);
        for record in &records {
            let obs = sample.fragment_observation(record, 60u8, variant);
            println!("{:?}", obs);
            assert_relative_eq!(obs.prob_ref.exp(), 0.0, epsilon=0.001);
            assert!(obs.prob_alt > obs.prob_ref);
        }
    }
}

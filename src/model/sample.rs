use std::str;
use std::collections::{HashMap, VecDeque, vec_deque};
use std::cmp;
use std::error::Error;
use std::f64::consts;
use std::cell::RefCell;
use std::io;
use std::f64;
use std::str::FromStr;

use csv;
use ordered_float::NotNaN;
use rgsl::randist::gaussian::{gaussian_pdf, ugaussian_P};
use rgsl::error::erfc;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::CigarStringView;
use bio::stats::{LogProb, Prob};

use model;
use model::Variant;
use model::evidence;
use model::evidence::{observation, Observation, Evidence};
use utils::Overlap;


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
    pub window: u32,
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
        if let Some(tid) = self.reader.header().tid(chrom) {
            let window_start = cmp::max(start as i32 - self.window as i32 - 1, 0) as u32;
            if self.inner.is_empty() || self.end().unwrap() < window_start || self.tid().unwrap() != tid as i32 {
                let end = self.reader.header().target_len(tid).unwrap();
                try!(self.reader.fetch(tid, window_start, end));
                debug!("Clearing ringbuffer");
                self.inner.clear();
            } else {
                // remove records too far left
                let to_remove = self.inner.iter().take_while(|rec| rec.pos() < window_start as i32).count();
                debug!("Removing {} records", to_remove);
                for _ in 0..to_remove {
                    self.inner.pop_front();
                }
            }

            // extend to the right
            loop {
                let mut record = bam::Record::new();
                if let Err(e) = self.reader.read(&mut record) {
                    if e.is_eof() {
                        break;
                    }
                    return Err(Box::new(e));
                }

                let pos = record.pos();
                if record.is_duplicate() || record.is_unmapped() || record.is_quality_check_failed() {
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

            debug!("New buffer length: {}", self.inner.len());

            Ok(())
        } else {
            Err(Box::new(RecordBufferError::UnknownSequence(str::from_utf8(chrom).unwrap().to_owned())))
        }
    }

    pub fn iter(&self) -> vec_deque::Iter<bam::Record> {
        self.inner.iter()
    }
}


struct Candidate<'a> {
    left: &'a bam::Record,
    right: Option<&'a bam::Record>
}


impl<'a> Candidate<'a> {
    fn new(record: &'a bam::Record) -> Self {
        Candidate {
            left: record,
            right: None
        }
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


/// A sequenced sample, e.g., a tumor or a normal sample.
pub struct Sample {
    record_buffer: RecordBuffer,
    use_fragment_evidence: bool,
    use_mapq: bool,
    adjust_mapq: bool,
    likelihood_model: model::likelihood::LatentVariableModel,
    pub(crate) indel_read_evidence: RefCell<evidence::reads::IndelEvidence>,
    pub(crate) indel_fragment_evidence: RefCell<evidence::fragments::IndelEvidence>
}


impl Sample {
    /// Create a new `Sample`.
    ///
    /// # Arguments
    ///
    /// * `bam` - BAM file with the aligned and deduplicated sequence reads.
    /// * `pileup_window` - Window around the variant that shall be searched for evidence (e.g. 5000).
    /// * `use_fragment_evidence` - Whether to use read pairs that are left and right of variant.
    /// * `use_secondary` - Whether to use secondary alignments.
    /// * `insert_size` - estimated insert size
    /// * `prior_model` - Prior assumptions about allele frequency spectrum of this sample.
    /// * `likelihood_model` - Latent variable model to calculate likelihoods of given observations.
    /// * `max_indel_overlap` - maximum number of bases a read may be aligned beyond the start or end of an indel in order to be considered as an observation
    /// * `indel_haplotype_window` - maximum number of considered bases around an indel breakpoint
    pub fn new(
        bam: bam::IndexedReader,
        pileup_window: u32,
        use_fragment_evidence: bool,
        // TODO remove this parameter, it will lead to wrong insert size estimations and is not necessary
        use_secondary: bool,
        // TODO remove this parameter, we should always use MAPQ
        use_mapq: bool,
        // TODO remove this parameter, it is not needed anymore
        adjust_mapq: bool,
        insert_size: InsertSize,
        likelihood_model: model::likelihood::LatentVariableModel,
        prob_insertion_artifact: Prob,
        prob_deletion_artifact: Prob,
        prob_insertion_extend_artifact: Prob,
        prob_deletion_extend_artifact: Prob,
        indel_haplotype_window: u32
    ) -> Self {
        Sample {
            record_buffer: RecordBuffer::new(bam, pileup_window, use_secondary),
            use_fragment_evidence: use_fragment_evidence,
            use_mapq: use_mapq,
            adjust_mapq: adjust_mapq,
            likelihood_model: likelihood_model,
            indel_read_evidence: RefCell::new(evidence::reads::IndelEvidence::new(
                LogProb::from(prob_insertion_artifact),
                LogProb::from(prob_deletion_artifact),
                LogProb::from(prob_insertion_extend_artifact),
                LogProb::from(prob_deletion_extend_artifact),
                indel_haplotype_window
            )),
            indel_fragment_evidence: RefCell::new(evidence::fragments::IndelEvidence::new(
                insert_size
            ))
        }
    }

    /// Return true if MAPQ appears to be reliable.
    /// Currently, this checks if AS > XS, i.e., the alignment score of the current position is
    /// better than for any alternative hit. If this is not the case, the read was most likely
    /// mapped to the current position because of its mate. Such placements can easily lead to
    /// false positives, especially in repetetive regions. Hence, we choose to rather ignore them.
    fn is_reliable_read(&self, record: &bam::Record) -> bool {
        if let Some(astag) = record.aux(b"AS") {
            if let Some(xstag) = record.aux(b"XS") {
                return astag.integer() > xstag.integer();
            }
        }

        true
    }

    /// Return likelihood model.
    pub fn likelihood_model(&self) -> model::likelihood::LatentVariableModel {
        self.likelihood_model
    }

    pub fn extract_common_observations(
        &mut self,
        chrom: &[u8],
        start: u32,
        variant: &Variant,
        common_obs: &mut observation::Common
    ) -> Result<(), Box<Error>> {
        let end = variant.end(start);
        // move window to the current variant
        self.record_buffer.fill(chrom, start, end)?;

        // Obtain common observations.
        common_obs.update(&self.record_buffer, start, variant)?;

        Ok(())
    }

    /// Extract observations for the given variant.
    pub fn extract_observations(
        &mut self,
        start: u32,
        variant: &Variant,
        chrom_name: &[u8],
        chrom_seq: &[u8],
        common_obs: &observation::Common
    ) -> Result<Vec<Observation>, Box<Error>> {
        let centerpoint = variant.centerpoint(start);

        let mut observations = Vec::new();
        let mut candidate_records = HashMap::new();

        match variant {
            //TODO: make &Variant::Ref add reads with position deleted if we want to check against indel alt alleles
            &Variant::SNV(_) | &Variant::None => {
                // iterate over records
                for record in self.record_buffer.iter() {
                    // TODO remove
                    if !self.is_reliable_read(record) {
                        continue;
                    }
                    let cigar = record.cigar();
                    let overlap = Overlap::new(
                        record, &cigar, start, variant, false
                    )?;

                    if overlap.is_enclosing() {
                        if let Some( obs ) = self.read_observation(
                            &record, &cigar, start, variant, chrom_name, chrom_seq, common_obs
                        )? {
                            observations.push( obs );
                        } else {
                             debug!("Did not add read to observations, SNV position deleted (Cigar op 'D') or skipped (Cigar op 'N').");
                        }
                    }
                }
            },
            &Variant::Insertion(_) | &Variant::Deletion(_) => {
                // iterate over records
                for record in self.record_buffer.iter() {

                    // We look at the whole fragment at once.

                    // We ensure fair sampling by checking if the whole fragment overlaps the
                    // centerpoint. Only taking the internal segment would not be fair,
                    // because then the second read of reference fragments tends to cross
                    // the centerpoint and the fragment would be discarded.
                    // The latter would not happen for alt (deletion) fragments, because the second
                    // read would map right of the variant in that case.

                    // We always choose the leftmost and the rightmost alignment, thereby also
                    // considering supplementary alignments.
                    if !candidate_records.contains_key(record.qname())
                    {
                        // this is the first (primary or supplementary alignment in the pair
                        candidate_records.insert(
                            record.qname().to_owned(),
                            Candidate::new(record)
                        );
                    } else if let Some(candidate) = candidate_records.get_mut(record.qname()) {
                        // this is either the last alignment or one in the middle
                        if (
                            candidate.left.is_first_in_template() && record.is_first_in_template()
                        ) && (
                            candidate.left.is_last_in_template() && record.is_last_in_template()
                        ) {
                            // ignore another partial alignment right of the first
                            continue;
                        }
                        // replace right record (we seek for the rightmost (partial) alignment)
                        candidate.right = Some(record);
                    }
                }
                for candidate in candidate_records.values() {
                    if let Some(right) = candidate.right {
                        // this is a pair
                        let start_pos = (candidate.left.pos() as u32).saturating_sub(
                            evidence::Clips::leading(&candidate.left.cigar()).soft()
                        );
                        if start_pos > centerpoint {
                            // ignore fragments that start beyond the centerpoint
                            continue;
                        }

                        let cigar = right.cigar();
                        let end_pos = cigar.end_pos()? as u32 +
                                      evidence::Clips::trailing(&cigar).soft();

                        if end_pos < centerpoint {
                            continue;
                        }
                        if let Some(obs) = self.fragment_observation(
                            candidate.left, right, start, variant, chrom_name, chrom_seq,
                            common_obs
                        )? {
                            observations.push(obs);
                        }
                    } else {
                        // this is a single alignment with unmapped mate or mate outside of the
                        // region of interest
                        let cigar = candidate.left.cigar();
                        let overlap = Overlap::new(
                            candidate.left, &cigar, start, variant, true
                        )?;
                        if !overlap.is_none() && candidate.left.is_mate_unmapped() {
                            if let Some(obs) = self.read_observation(
                                    candidate.left, &cigar, start, variant, chrom_name,
                                    chrom_seq, common_obs
                            )? {
                                observations.push(obs);
                            }
                        }
                    }
                }
            }
        }

        if !observations.is_empty() {
            // We scale all probabilities by the maximum value. This is just an unbiased scaling
            // that does not affect the final certainties (because of Bayes' theorem application
            // in the end). However, we avoid numerical issues (e.g., during integration).
            let max_prob = LogProb(*observations.iter().map(|obs| {
                cmp::max(NotNaN::from(obs.prob_ref), NotNaN::from(obs.prob_alt))
            }).max().unwrap());
            if max_prob != LogProb::ln_zero() {
                // only scale if the maximum probability is not zero
                for obs in observations.iter_mut() {
                    obs.prob_ref = obs.prob_ref - max_prob;
                    obs.prob_alt = obs.prob_alt - max_prob;
                    assert!(obs.prob_ref.is_valid());
                    assert!(obs.prob_alt.is_valid());
                }
            }
        }
        Ok(observations)
    }

    /// Obtain mapping quality. If adjust_mapq is set, the mapping quality is adjusted with the
    /// certainty of the aligment score compared to the best alternative alignment score.
    fn prob_mapping(
        &self,
        record: &bam::Record,
        cigar: &bam::record::CigarStringView,
        chrom_name: &[u8],
        chrom_seq: &[u8]
    ) -> Result<LogProb, Box<Error>> {
        if self.use_mapq {
            if self.adjust_mapq {
                Ok(evidence::reads::prob_mapping_adjusted(
                    record,
                    cigar,
                    chrom_name,
                    chrom_seq
                )?)
            } else {
                Ok(evidence::reads::prob_mapping(record))
            }
        } else {
            Ok(LogProb::ln_one())
        }
    }

    /// extract within-read evidence for reads covering an indel or SNV of interest
    fn read_observation(
        &self,
        record: &bam::Record,
        cigar: &CigarStringView,
        start: u32,
        variant: &Variant,
        chrom_name: &[u8],
        chrom_seq: &[u8],
        common_obs: &observation::Common
    ) -> Result<Option<Observation>, Box<Error>> {
        let probs = match variant {
            &Variant::Deletion(_) | &Variant::Insertion(_) => {
                Some( self.indel_read_evidence.borrow_mut()
                                        .prob(record, cigar, start, variant, chrom_seq)? )
            },
            &Variant::SNV(_) => {
                evidence::reads::prob_snv(record, &cigar, start, variant, chrom_seq)?
            },
            &Variant::None => {
                evidence::reads::prob_none(record, &cigar, start, variant, chrom_seq)?
            }
        };

        if let Some( (prob_ref, prob_alt) ) = probs {
            let prob_mapping = self.prob_mapping(record, cigar, chrom_name, chrom_seq)?;

            let prob_sample_alt = self.indel_read_evidence.borrow().prob_sample_alt(
                record.seq().len() as u32,
                variant
            );
            Ok( Some (
                Observation::new(
                    prob_mapping,
                    prob_alt,
                    prob_ref,
                    prob_sample_alt.joint_prob(common_obs),
                    Evidence::alignment(cigar, record)
                )
            ))
        } else {
            Ok( None )
        }
    }

    /// Extract insert size information for fragments (e.g. read pairs) spanning an indel of interest
    /// Here we calculate the product of insert size based and alignment based probabilities.
    /// This has the benefit that the calculation automatically checks for consistence between
    /// insert size and overlapping alignmnments.
    /// This sports the following desirable behavior:
    ///
    /// * If there is no clear evidence from either the insert size or the alignment, the factors
    ///   simply scale because the probabilities of the corresponding type of evidence will be equal.
    /// * If reads and fragments agree, 1 stays 1 and 0 stays 0.
    /// * If reads and fragments disagree (the promising part!), the other type of evidence will
    ///   scale potential false positive probabilities down.
    /// * Since there is only one observation per fragment, there is no double counting when
    ///   estimating allele frequencies. Before, we had one observation for an overlapping read
    ///   and potentially another observation for the corresponding fragment.
    fn fragment_observation(
        &self,
        left_record: &bam::Record,
        right_record: &bam::Record,
        start: u32,
        variant: &Variant,
        chrom_name: &[u8],
        chrom_seq: &[u8],
        common_obs: &observation::Common
    ) -> Result<Option<Observation>, Box<Error>> {

        let prob_read = |
            record: &bam::Record, cigar: &CigarStringView
        | -> Result<(LogProb, LogProb), Box<Error>> {
            // Calculate read evidence.
            // We also calculate it in case of no overlap. Otherwise, there would be a bias due to
            // non-overlapping fragments having higher likelihoods.
            Ok(self.indel_read_evidence.borrow_mut()
                                       .prob(record, cigar, start, variant, chrom_seq)?)
        };

        let left_cigar = left_record.cigar();
        let right_cigar = right_record.cigar();
        let left_overlap = Overlap::new(
            left_record, &left_cigar, start, variant, true
        )?;
        let right_overlap = Overlap::new(
            right_record, &right_cigar, start, variant, true
        )?;

        if !self.use_fragment_evidence && left_overlap.is_none() && right_overlap.is_none() {
            return Ok(None);
        }

        let (p_ref_left, p_alt_left) = prob_read(left_record, &left_cigar)?;
        let (p_ref_right, p_alt_right) = prob_read(right_record, &right_cigar)?;

        let left_read_len = left_record.seq().len() as u32;
        let right_read_len = right_record.seq().len() as u32;

        let insert_size = evidence::fragments::estimate_insert_size(left_record, right_record)?;
        let (p_ref_isize, p_alt_isize) = match variant {
            &Variant::Deletion(_) if self.use_fragment_evidence => {
                self.indel_fragment_evidence.borrow().prob(insert_size, variant)?
            },
            _ => {
                // Ignore isize for insertions. The reason is that we cannot reliably determine if a
                // fragment encloses the insertion properly (with overlaps at the inner read ends).
                // Hence, the probabilities cannot be calculated. Further, we have a lot of fragments
                // that overlap insertions at the left or right side, and those are also helpful.
                (LogProb::ln_one(), LogProb::ln_one())
            }
        };

        let prob_sample_alt = self.indel_fragment_evidence.borrow().prob_sample_alt(
            left_read_len,
            right_read_len,
            variant
        );

        assert!(p_alt_isize.is_valid());
        assert!(p_ref_isize.is_valid());
        assert!(p_alt_left.is_valid());
        assert!(p_alt_right.is_valid());
        assert!(p_ref_left.is_valid());
        assert!(p_ref_right.is_valid());

        let obs = Observation::new(
            self.prob_mapping(left_record, &left_cigar, chrom_name, chrom_seq)? +
            self.prob_mapping(right_record, &right_cigar, chrom_name, chrom_seq)?,
            p_alt_isize + p_alt_left + p_alt_right,
            p_ref_isize + p_ref_left + p_ref_right,
            prob_sample_alt.joint_prob(common_obs),
            Evidence::insert_size(
                insert_size as u32,
                &left_record.cigar(),
                &right_record.cigar(),
                left_record,
                right_record,
                p_ref_left,
                p_alt_left,
                p_ref_right,
                p_alt_right,
                p_ref_isize,
                p_alt_isize
            )
        );
        assert!(obs.prob_alt.is_valid());
        assert!(obs.prob_ref.is_valid());

        Ok(Some(obs))
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
    extern crate env_logger;

    use super::*;
    use model;
    use likelihood;
    use constants;

    use std::str;
    use std::fs;
    use itertools::Itertools;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use bio::stats::{LogProb, PHREDProb, Prob};
    use bio::io::fasta;
    use model::tests::common_observation;


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

    fn setup_sample(isize_mean: f64) -> Sample {
        Sample::new(
            bam::IndexedReader::from_path(&"tests/indels.bam").unwrap(),
            2500,
            true,
            true,
            true,
            false,
            InsertSize { mean: isize_mean, sd: 20.0 },
            likelihood::LatentVariableModel::new(1.0),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            10
        )
    }

    #[test]
    #[ignore]
    fn test_read_observation_indel() {
        let variant = model::Variant::Insertion(b"GCATCCTGCG".to_vec());
        // insertion starts at 546 and has length 10
        let varpos = 546;

        let sample = setup_sample(150.0);
        let mut bam = bam::Reader::from_path(&"tests/indels.bam").unwrap();
        let records = bam.records().collect_vec();

        let ref_seq = ref_seq();

        let true_alt_probs = [-0.09, -0.02, -73.09, -16.95, -73.09];

        let common = common_observation();

        for (record, true_alt_prob) in records.into_iter().zip(true_alt_probs.into_iter()) {
            let record = record.unwrap();
            let cigar = record.cigar();
            if let Some( obs ) = sample.read_observation(&record, &cigar, varpos, &variant, b"17", &ref_seq, &common).unwrap() {
                println!("{:?}", obs);
                assert_relative_eq!(*obs.prob_alt, *true_alt_prob, epsilon=0.01);
                assert_relative_eq!(*obs.prob_mapping, *(LogProb::from(PHREDProb(60.0)).ln_one_minus_exp()));
            } else {
                panic!("read_observation() test for indels failed; it returned 'None'.")
            }
        }
    }

    #[test]
    fn test_record_buffer() {
        let bam = bam::IndexedReader::from_path(&"tests/indels.bam").unwrap();
        let mut buffer = RecordBuffer::new(bam, 10, true);

        buffer.fill(b"17", 10, 20).unwrap();
        buffer.fill(b"17", 478, 500).unwrap();
        buffer.fill(b"17", 1000, 1700).unwrap();
        // TODO add assertions
    }

    fn ref_seq() -> Vec<u8> {
        let mut fa = fasta::Reader::from_file(&"tests/chr17.prefix.fa").unwrap();
        let mut chr17 = fasta::Record::new();
        fa.read(&mut chr17).unwrap();

        chr17.seq().to_owned()
    }

    #[test]
    fn test_prob_read_indel() {
        let _ = env_logger::init();

        let mut bam = bam::Reader::from_path(&"tests/indels+clips.bam").unwrap();
        let records = bam.records().map(|rec| rec.unwrap()).collect_vec();
        let ref_seq = ref_seq();
        let sample = setup_sample(312.0);

        // truth
        let probs_alt = [-0.09, -16.95, -73.09, -0.022, -0.011, -0.03];
        let probs_ref = [-151.13, -163.03, -0.01, -67.75, -67.74, -67.76];

        // variant (obtained via bcftools)
        let start = 546;
        let variant = model::Variant::Insertion(b"GCATCCTGCG".to_vec());
        for (i, rec) in records.iter().enumerate() {
            println!("{}", str::from_utf8(rec.qname()).unwrap());
            let (prob_ref, prob_alt) = sample.indel_read_evidence.borrow_mut().prob(rec, &rec.cigar(), start, &variant, &ref_seq).unwrap();
            println!("Pr(ref)={} Pr(alt)={}", *prob_ref, *prob_alt);
            println!("{:?}", rec.cigar());
            assert_relative_eq!(*prob_ref, probs_ref[i], epsilon=0.1);
            assert_relative_eq!(*prob_alt, probs_alt[i], epsilon=0.1);
        }
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

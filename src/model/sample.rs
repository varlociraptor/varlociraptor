use std::cell::{RefCell, RefMut};
use std::cmp;
use std::collections::{vec_deque, BTreeMap, VecDeque};
use std::error::Error;
use std::f64;
use std::f64::consts;
use std::str;

use bio::stats::{LogProb, Prob};
use rand::distributions;
use rand::distributions::IndependentSample;
use rand::{SeedableRng, StdRng};
use rgsl::error::erfc;
use rgsl::randist::gaussian::{gaussian_pdf, ugaussian_P};
use rust_htslib::bam;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::bam::Read;

use crate::estimation::alignment_properties;
use crate::model::evidence;
use crate::model::evidence::reads::AbstractReadEvidence;
use crate::model::evidence::{Evidence, Observation};
use crate::model::{Variant, VariantType};
use crate::utils::{is_repeat_variant, max_prob, Overlap};

quick_error! {
    #[derive(Debug)]
    pub enum RecordBufferError {
        UnknownSequence(chrom: String) {
            description("unknown sequence")
            display("sequence {} cannot be found in BAM", chrom)
        }
    }
}

pub type Pileup = Vec<Observation>;

/// Ringbuffer of BAM records. This data structure ensures that no bam record is read twice while
/// extracting observations for given variants.
pub struct RecordBuffer {
    reader: bam::IndexedReader,
    inner: VecDeque<bam::Record>,
    pub window: u32,
    use_secondary: bool,
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
            use_secondary: use_secondary,
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
            if self.inner.is_empty()
                || self.end().unwrap() < window_start
                || self.tid().unwrap() != tid as i32
            {
                let end = self.reader.header().target_len(tid).unwrap();
                self.reader.fetch(tid, window_start, end)?;
                debug!("Clearing ringbuffer");
                self.inner.clear();
            } else {
                // remove records too far left
                let to_remove = self
                    .inner
                    .iter()
                    .take_while(|rec| rec.pos() < window_start as i32)
                    .count();
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
                if record.is_duplicate() || record.is_unmapped() || record.is_quality_check_failed()
                {
                    continue;
                }
                if !self.use_secondary && record.is_secondary() {
                    continue;
                }
                // unpack cigar string
                record.cache_cigar();
                self.inner.push_back(record);
                if pos > end as i32 + self.window as i32 {
                    break;
                }
            }

            debug!("New buffer length: {}", self.inner.len());

            Ok(())
        } else {
            Err(Box::new(RecordBufferError::UnknownSequence(
                str::from_utf8(chrom).unwrap().to_owned(),
            )))
        }
    }

    pub fn iter(&self) -> vec_deque::Iter<bam::Record> {
        self.inner.iter()
    }
}

struct Candidate<'a> {
    left: &'a bam::Record,
    right: Option<&'a bam::Record>,
}

impl<'a> Candidate<'a> {
    fn new(record: &'a bam::Record) -> Self {
        Candidate {
            left: record,
            right: None,
        }
    }
}

pub enum SubsampleCandidates {
    Necessary {
        rng: StdRng,
        prob: f64,
        prob_range: distributions::Range<f64>,
    },
    None,
}

impl SubsampleCandidates {
    pub fn new(max_depth: usize, depth: usize) -> Self {
        if depth > max_depth {
            SubsampleCandidates::Necessary {
                rng: StdRng::from_seed(&[48074578]),
                prob: max_depth as f64 / depth as f64,
                prob_range: distributions::Range::new(0.0, 1.0),
            }
        } else {
            SubsampleCandidates::None
        }
    }

    pub fn keep(&mut self) -> bool {
        match self {
            SubsampleCandidates::Necessary {
                rng,
                prob,
                prob_range,
            } => prob_range.ind_sample(rng) <= *prob,
            SubsampleCandidates::None => true,
        }
    }
}

/// A sequenced sample, e.g., a tumor or a normal sample.
pub struct Sample {
    record_buffer: RecordBuffer,
    use_fragment_evidence: bool,
    alignment_properties: alignment_properties::AlignmentProperties,
    pub(crate) indel_read_evidence: RefCell<evidence::reads::IndelEvidence>,
    pub(crate) indel_fragment_evidence: RefCell<evidence::fragments::IndelEvidence>,
    pub(crate) snv_read_evidence: RefCell<evidence::reads::SNVEvidence>,
    pub(crate) none_read_evidence: RefCell<evidence::reads::NoneEvidence>,
    max_depth: usize,
    omit_repeat_regions: Vec<VariantType>,
}

impl Sample {
    /// Create a new `Sample`.
    ///
    /// # Arguments
    ///
    /// * `bam` - BAM file with the aligned and deduplicated sequence reads.
    /// * `use_fragment_evidence` - Whether to use read pairs that are left and right of variant.
    /// * `use_secondary` - Whether to use secondary alignments.
    /// * `insert_size` - estimated insert size
    /// * `prior_model` - Prior assumptions about allele frequency spectrum of this sample.
    /// * `max_indel_overlap` - maximum number of bases a read may be aligned beyond the start or end of an indel in order to be considered as an observation
    /// * `indel_haplotype_window` - maximum number of considered bases around an indel breakpoint
    pub fn new(
        bam: bam::IndexedReader,
        use_fragment_evidence: bool,
        alignment_properties: alignment_properties::AlignmentProperties,
        prob_insertion_artifact: Prob,
        prob_deletion_artifact: Prob,
        prob_insertion_extend_artifact: Prob,
        prob_deletion_extend_artifact: Prob,
        indel_haplotype_window: u32,
        max_depth: usize,
        omit_repeat_regions: &[VariantType],
    ) -> Self {
        let pileup_window = (alignment_properties.insert_size().mean
            + alignment_properties.insert_size().sd * 6.0) as u32;
        info!(
            "Using window of {} bases on each side of variant.",
            pileup_window
        );

        Sample {
            record_buffer: RecordBuffer::new(bam, pileup_window, false),
            use_fragment_evidence: use_fragment_evidence,
            alignment_properties: alignment_properties,
            indel_read_evidence: RefCell::new(evidence::reads::IndelEvidence::new(
                LogProb::from(prob_insertion_artifact),
                LogProb::from(prob_deletion_artifact),
                LogProb::from(prob_insertion_extend_artifact),
                LogProb::from(prob_deletion_extend_artifact),
                indel_haplotype_window,
            )),
            snv_read_evidence: RefCell::new(evidence::reads::SNVEvidence::new()),
            indel_fragment_evidence: RefCell::new(evidence::fragments::IndelEvidence::new()),
            none_read_evidence: RefCell::new(evidence::reads::NoneEvidence::new()),
            max_depth: max_depth,
            omit_repeat_regions: omit_repeat_regions.to_vec(),
        }
    }

    /// Extract observations for the given variant.
    pub fn extract_observations(
        &mut self,
        start: u32,
        variant: &Variant,
        chrom: &[u8],
        chrom_seq: &[u8],
    ) -> Result<Pileup, Box<Error>> {
        let centerpoint = variant.centerpoint(start);

        for vartype in &self.omit_repeat_regions {
            if variant.is_type(vartype) && is_repeat_variant(start, variant, chrom_seq) {
                // Do not return evidence, in order to mark variant in output as unclear.
                return Ok(Vec::new());
            }
        }

        self.record_buffer.fill(chrom, start, variant.end(start))?;

        let mut observations = Vec::new();

        match variant {
            //TODO: make &Variant::None add reads with position deleted if we want to check against indel alt alleles
            &Variant::SNV(_) | &Variant::None => {
                let mut candidate_records = Vec::new();
                // iterate over records
                for record in self.record_buffer.iter() {
                    if record.pos() as u32 > start {
                        // the read cannot overlap the variant
                        continue;
                    }

                    let cigar = record.cigar_cached().unwrap();
                    let overlap = Overlap::new(record, cigar, start, variant, false)?;

                    if overlap.is_enclosing() {
                        candidate_records.push(record);
                    }
                }
                let mut subsample_candidates =
                    SubsampleCandidates::new(self.max_depth, candidate_records.len());

                for record in candidate_records {
                    if subsample_candidates.keep() {
                        if let Some(obs) = self.read_observation(
                            &record,
                            record.cigar_cached().unwrap(),
                            start,
                            variant,
                            chrom_seq,
                        )? {
                            observations.push(obs);
                        } else {
                            debug!("Did not add read to observations, SNV position deleted (Cigar op 'D') or skipped (Cigar op 'N').");
                        }
                    }
                }
            }
            &Variant::Insertion(_) | &Variant::Deletion(_) => {
                // We cannot use a hash function here because candidates have to be considered
                // in a deterministic order. Otherwise, subsampling high-depth regions will result
                // in slightly different probabilities each time.
                let mut candidate_records = BTreeMap::new();

                // iterate over records
                for record in self.record_buffer.iter() {
                    // First, we check whether the record contains an indel in the cigar.
                    // We store the maximum indel size to update the global estimates, in case
                    // it is larger in this region.
                    self.alignment_properties.update_max_cigar_ops_len(record);

                    // We look at the whole fragment at once.

                    // We ensure fair sampling by checking if the whole fragment overlaps the
                    // centerpoint. Only taking the internal segment would not be fair,
                    // because then the second read of reference fragments tends to cross
                    // the centerpoint and the fragment would be discarded.
                    // The latter would not happen for alt (deletion) fragments, because the second
                    // read would map right of the variant in that case.

                    // We always choose the leftmost and the rightmost alignment, thereby also
                    // considering supplementary alignments.
                    if !candidate_records.contains_key(record.qname()) {
                        // this is the first (primary or supplementary alignment in the pair
                        candidate_records.insert(record.qname().to_owned(), Candidate::new(record));
                    } else if let Some(candidate) = candidate_records.get_mut(record.qname()) {
                        // this is either the last alignment or one in the middle
                        if (candidate.left.is_first_in_template() && record.is_first_in_template())
                            && (candidate.left.is_last_in_template()
                                && record.is_last_in_template())
                        {
                            // ignore another partial alignment right of the first
                            continue;
                        }
                        // replace right record (we seek for the rightmost (partial) alignment)
                        candidate.right = Some(record);
                    }
                }

                let mut candidate_fragments = Vec::new();
                let mut candidate_reads = Vec::new();
                for candidate in candidate_records.values() {
                    if let Some(right) = candidate.right {
                        // this is a pair
                        let start_pos = (candidate.left.pos() as u32).saturating_sub(
                            evidence::Clips::leading(candidate.left.cigar_cached().unwrap()).soft(),
                        );
                        if start_pos > centerpoint {
                            // ignore fragments that start beyond the centerpoint
                            continue;
                        }

                        let cigar = right.cigar_cached().unwrap();
                        let end_pos =
                            cigar.end_pos()? as u32 + evidence::Clips::trailing(cigar).soft();

                        if end_pos < centerpoint {
                            continue;
                        }

                        let left_cigar = candidate.left.cigar_cached().unwrap();
                        let right_cigar = right.cigar_cached().unwrap();
                        let left_overlap =
                            Overlap::new(candidate.left, left_cigar, start, variant, true)?;
                        let right_overlap = Overlap::new(right, right_cigar, start, variant, true)?;

                        if left_overlap.is_none() && right_overlap.is_none() {
                            // Skip fragment if none of the reads overlaps the variant.
                            // This increases robustness, because insert size is never considered alone.
                            continue;
                        }

                        candidate_fragments.push(candidate);
                    } else {
                        // this is a single alignment with unmapped mate or mate outside of the
                        // region of interest
                        let cigar = candidate.left.cigar_cached().unwrap();
                        let overlap = Overlap::new(candidate.left, cigar, start, variant, true)?;
                        if !overlap.is_none() && candidate.left.is_mate_unmapped() {
                            candidate_reads.push(candidate);
                        }
                    }
                }

                let mut subsample_candidates = SubsampleCandidates::new(
                    self.max_depth,
                    candidate_fragments.len() + candidate_reads.len(),
                );

                for candidate in candidate_fragments {
                    if !subsample_candidates.keep() {
                        continue;
                    }

                    if let Some(obs) = self.fragment_observation(
                        candidate.left,
                        candidate.right.unwrap(),
                        start,
                        variant,
                        chrom_seq,
                    )? {
                        observations.push(obs);
                    }
                }

                for candidate in candidate_reads {
                    if !subsample_candidates.keep() {
                        continue;
                    }
                    if let Some(obs) = self.read_observation(
                        candidate.left,
                        candidate.left.cigar_cached().unwrap(),
                        start,
                        variant,
                        chrom_seq,
                    )? {
                        observations.push(obs);
                    }
                }

                //self.refine_prob_mapping(&mut observations);
            }
        }

        Ok(observations)
    }

    /// extract within-read evidence for reads covering an indel or SNV of interest
    fn read_observation(
        &self,
        record: &bam::Record,
        cigar: &CigarStringView,
        start: u32,
        variant: &Variant,
        chrom_seq: &[u8],
    ) -> Result<Option<Observation>, Box<Error>> {
        let mut evidence: RefMut<evidence::reads::AbstractReadEvidence> = match variant {
            &Variant::Deletion(_) | &Variant::Insertion(_) => self.indel_read_evidence.borrow_mut(),
            &Variant::SNV(_) => self.snv_read_evidence.borrow_mut(),
            &Variant::None => self.none_read_evidence.borrow_mut(),
        };

        if let Some((prob_ref, prob_alt)) =
            evidence.prob(record, cigar, start, variant, chrom_seq)?
        {
            let (prob_mapping, prob_mismapping) = evidence.prob_mapping_mismapping(record);

            // This is an estimate of the allele likelihood at the true location in case the read is
            // mismapped.
            let prob_missed_allele = max_prob(prob_ref, prob_alt);

            let prob_sample_alt = evidence.prob_sample_alt(
                record.seq().len() as u32,
                variant,
                &self.alignment_properties,
            );
            Ok(Some(Observation::new(
                prob_mapping,
                prob_mismapping,
                prob_alt,
                prob_ref,
                prob_missed_allele,
                prob_sample_alt,
                Evidence::alignment(cigar, record),
            )))
        } else {
            Ok(None)
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
        chrom_seq: &[u8],
    ) -> Result<Option<Observation>, Box<Error>> {
        let prob_read = |record: &bam::Record,
                         cigar: &CigarStringView|
         -> Result<(LogProb, LogProb), Box<Error>> {
            // Calculate read evidence.
            // We also calculate it in case of no overlap. Otherwise, there would be a bias due to
            // non-overlapping fragments having higher likelihoods.
            Ok(self
                .indel_read_evidence
                .borrow_mut()
                .prob(record, cigar, start, variant, chrom_seq)?
                .unwrap())
        };

        let left_cigar = left_record.cigar_cached().unwrap();
        let right_cigar = right_record.cigar_cached().unwrap();

        let (p_ref_left, p_alt_left) = prob_read(left_record, left_cigar)?;
        let (p_ref_right, p_alt_right) = prob_read(right_record, right_cigar)?;

        // This is an estimate of the allele likelihood at the true location in case the read is
        // mismapped.
        let p_missed_left = max_prob(p_ref_left, p_alt_left);
        let p_missed_right = max_prob(p_ref_right, p_alt_right);

        let left_read_len = left_record.seq().len() as u32;
        let right_read_len = right_record.seq().len() as u32;

        let insert_size = evidence::fragments::estimate_insert_size(left_record, right_record)?;
        let (p_ref_isize, p_alt_isize) = match variant {
            &Variant::Deletion(_) if self.use_fragment_evidence => self
                .indel_fragment_evidence
                .borrow()
                .prob(insert_size, variant, &self.alignment_properties)?,
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
            variant,
            &self.alignment_properties,
        );

        assert!(p_alt_isize.is_valid());
        assert!(p_ref_isize.is_valid());
        assert!(p_alt_left.is_valid());
        assert!(p_alt_right.is_valid());
        assert!(p_ref_left.is_valid());
        assert!(p_ref_right.is_valid());

        let (_, prob_mismapping_left) = self
            .indel_read_evidence
            .borrow()
            .prob_mapping_mismapping(left_record);
        let (_, prob_mismapping_right) = self
            .indel_read_evidence
            .borrow()
            .prob_mapping_mismapping(right_record);
        let prob_mismapping = prob_mismapping_left + prob_mismapping_right;
        let obs = Observation::new(
            prob_mismapping.ln_one_minus_exp(),
            prob_mismapping,
            p_alt_isize + p_alt_left + p_alt_right,
            p_ref_isize + p_ref_left + p_ref_right,
            p_missed_left + p_missed_right,
            prob_sample_alt,
            Evidence::insert_size(
                insert_size as u32,
                left_record.cigar_cached().unwrap(),
                right_record.cigar_cached().unwrap(),
                left_record,
                right_record,
                p_ref_left,
                p_alt_left,
                p_ref_right,
                p_alt_right,
                p_ref_isize,
                p_alt_isize,
            ),
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
        (ugaussian_P((value + 0.5 - mean) / sd) - ugaussian_P((value - 0.5 - mean) / sd)).ln(), // - ugaussian_P(-mean / sd).ln()
    )
}

/// Continuous normal density (obsolete).
pub fn isize_density(value: f64, mean: f64, sd: f64) -> LogProb {
    LogProb(gaussian_pdf(value - mean, sd).ln())
}

/// Manual normal density (obsolete, we can use GSL (see above)).
pub fn isize_density_louis(value: f64, mean: f64, sd: f64) -> LogProb {
    let mut p = 0.5 / (1.0 - 0.5 * erfc((mean + 0.5) / sd * consts::FRAC_1_SQRT_2));
    p *= erfc((-value - 0.5 + mean) / sd * consts::FRAC_1_SQRT_2)
        - erfc((-value + 0.5 + mean) / sd * consts::FRAC_1_SQRT_2);

    LogProb(p.ln())
}

pub fn isize_mixture_density_louis(value: f64, d: f64, mean: f64, sd: f64, rate: f64) -> LogProb {
    let p = 0.5
        / (rate * (1.0 - 0.5 * erfc((mean + 0.5) / sd * consts::FRAC_1_SQRT_2))
            + (1.0 - rate) * (1.0 - 0.5 * erfc((mean + d + 0.5) / sd * consts::FRAC_1_SQRT_2)));
    LogProb(
        (p * (rate
            * (erfc((-value - 0.5 + mean) / sd * consts::FRAC_1_SQRT_2)
                - erfc((-value + 0.5 + mean) / sd * consts::FRAC_1_SQRT_2))
            + (1.0 - rate)
                * (erfc((-value - 0.5 + mean + d) / sd * consts::FRAC_1_SQRT_2)
                    - erfc((-value + 0.5 + mean + d) / sd * consts::FRAC_1_SQRT_2))))
        .ln(),
    )
}

#[cfg(test)]
mod tests {
    extern crate env_logger;

    use super::*;
    use crate::constants;
    use crate::likelihood;
    use crate::model;

    use bio::io::fasta::{self, FastaRead};
    use bio::stats::{LogProb, PHREDProb, Prob};
    use crate::estimation::alignment_properties::{AlignmentProperties, InsertSize};
    use itertools::Itertools;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::str;

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
            rate.ln_one_minus_exp() + isize_pmf(212.0, 312.0 - 100.0, 15.0)
        )
        .ln_add_exp(
            // case: no indel, false positive call
            rate + isize_pmf(212.0, 312.0, 15.0),
        );

        println!("{} {}", d_mix.exp(), p_alt.exp());
    }

    fn setup_sample(isize_mean: f64) -> Sample {
        Sample::new(
            bam::IndexedReader::from_path(&"tests/indels.bam").unwrap(),
            true,
            AlignmentProperties::default(InsertSize {
                mean: isize_mean,
                sd: 20.0,
            }),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            10,
            500,
            &[],
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

        for (record, true_alt_prob) in records.into_iter().zip(true_alt_probs.into_iter()) {
            let mut record = record.unwrap();
            record.cache_cigar();
            let cigar = record.cigar_cached().unwrap();
            if let Some(obs) = sample
                .read_observation(&record, cigar, varpos, &variant, &ref_seq)
                .unwrap()
            {
                println!("{:?}", obs);
                assert_relative_eq!(*obs.prob_alt, *true_alt_prob, epsilon = 0.01);
                assert_relative_eq!(
                    *obs.prob_mapping,
                    *(LogProb::from(PHREDProb(60.0)).ln_one_minus_exp())
                );
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

    // TODO re-enable once framework has stabilized
    #[test]
    #[ignore]
    fn test_prob_read_indel() {
        let _ = env_logger::init();

        let mut bam = bam::Reader::from_path(&"tests/indels+clips.bam").unwrap();
        let records = bam.records().map(|rec| rec.unwrap()).collect_vec();
        let ref_seq = ref_seq();
        let sample = setup_sample(312.0);

        // truth
        let probs_alt = [-0.09, -16.95, -73.09, -0.022, -0.011, -0.03];
        let probs_ref = [-150.50, -163.03, -0.01, -67.75, -67.74, -67.76];

        // variant (obtained via bcftools)
        let start = 546;
        let variant = model::Variant::Insertion(b"GCATCCTGCG".to_vec());
        for (i, mut rec) in records.into_iter().enumerate() {
            rec.cache_cigar();
            println!("{}", str::from_utf8(rec.qname()).unwrap());
            let (prob_ref, prob_alt) = sample
                .indel_read_evidence
                .borrow_mut()
                .prob(&rec, rec.cigar_cached().unwrap(), start, &variant, &ref_seq)
                .unwrap()
                .unwrap();
            println!("Pr(ref)={} Pr(alt)={}", *prob_ref, *prob_alt);
            println!("{:?}", rec.cigar_cached());
            assert_relative_eq!(*prob_ref, probs_ref[i], epsilon = 0.1);
            assert_relative_eq!(*prob_alt, probs_alt[i], epsilon = 0.1);
        }
    }
}

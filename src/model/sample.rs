// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::{RefCell, RefMut};
use std::collections::BTreeMap;
use std::error::Error;
use std::f64;
use std::path::Path;
use std::str;

use bio::stats::{LogProb, Prob};
use bio_types::strand::Strand;
use derive_builder::Builder;
use rand::distributions;
use rand::distributions::Distribution;
use rand::{rngs::StdRng, SeedableRng};
use rust_htslib::bam;
use rust_htslib::bam::record::CigarStringView;

use crate::estimation::alignment_properties;
use crate::model::evidence;
use crate::model::evidence::reads::AbstractReadEvidence;
use crate::model::evidence::{observation::ObservationBuilder, Observation};
use crate::model::{Variant, VariantType};
use crate::utils::{is_repeat_variant, Overlap};

/// Strand combination for read pairs as given by the sequencing protocol.
#[derive(
    Display, Debug, Clone, Copy, Serialize, Deserialize, EnumString, EnumIter, IntoStaticStr,
)]
pub enum ProtocolStrandedness {
    #[strum(serialize = "opposite")]
    Opposite,
    #[strum(serialize = "same")]
    Same,
}

impl Default for ProtocolStrandedness {
    fn default() -> Self {
        ProtocolStrandedness::Opposite
    }
}

pub type Pileup = Vec<Observation>;

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
        prob_range: distributions::Uniform<f64>,
    },
    None,
}

impl SubsampleCandidates {
    pub fn new(max_depth: usize, depth: usize) -> Self {
        if depth > max_depth {
            SubsampleCandidates::Necessary {
                rng: StdRng::seed_from_u64(48074578),
                prob: max_depth as f64 / depth as f64,
                prob_range: distributions::Uniform::new(0.0, 1.0),
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
            } => prob_range.sample(rng) <= *prob,
            SubsampleCandidates::None => true,
        }
    }
}

pub fn estimate_alignment_properties<P: AsRef<Path>>(
    path: P,
) -> Result<alignment_properties::AlignmentProperties, Box<dyn Error>> {
    let mut bam = bam::Reader::from_path(path)?;
    Ok(alignment_properties::AlignmentProperties::estimate(
        &mut bam,
    )?)
}

/// A sequenced sample, e.g., a tumor or a normal sample.
#[derive(Builder, Debug)]
#[builder(pattern = "owned")]
pub struct Sample {
    #[builder(private)]
    record_buffer: bam::buffer::RecordBuffer,
    #[builder(default = "true")]
    use_fragment_evidence: bool,
    #[builder(private)]
    alignment_properties: alignment_properties::AlignmentProperties,
    #[builder(private)]
    pub(crate) indel_read_evidence: RefCell<evidence::reads::IndelEvidence>,
    #[builder(private)]
    pub(crate) indel_fragment_evidence: RefCell<evidence::fragments::IndelEvidence>,
    #[builder(private)]
    pub(crate) snv_read_evidence: RefCell<evidence::reads::SNVEvidence>,
    #[builder(private)]
    pub(crate) mnv_read_evidence: RefCell<evidence::reads::MNVEvidence>,
    #[builder(private)]
    pub(crate) none_read_evidence: RefCell<evidence::reads::NoneEvidence>,
    #[builder(default = "200")]
    max_depth: usize,
    #[builder(default = "Vec::new()")]
    omit_repeat_regions: Vec<VariantType>,
    #[builder(private)]
    buffer_window: u32,
    protocol_strandedness: ProtocolStrandedness,
}

impl SampleBuilder {
    /// Register alignment information.
    ///
    /// # Arguments
    /// * `bam` - BAM file with the aligned and deduplicated sequence reads.
    pub fn alignments(
        self,
        bam: bam::IndexedReader,
        alignment_properties: alignment_properties::AlignmentProperties,
    ) -> Self {
        let pileup_window = (alignment_properties.insert_size().mean
            + alignment_properties.insert_size().sd * 6.0) as u32;
        self.alignment_properties(alignment_properties)
            .record_buffer(bam::RecordBuffer::new(bam, true))
            .buffer_window(pileup_window)
    }

    /// Register error probabilities and window to check around indels.
    pub fn error_probs(
        self,
        prob_insertion_artifact: Prob,
        prob_deletion_artifact: Prob,
        prob_insertion_extend_artifact: Prob,
        prob_deletion_extend_artifact: Prob,
        indel_haplotype_window: u32,
    ) -> Self {
        self.indel_read_evidence(RefCell::new(evidence::reads::IndelEvidence::new(
            LogProb::from(prob_insertion_artifact),
            LogProb::from(prob_deletion_artifact),
            LogProb::from(prob_insertion_extend_artifact),
            LogProb::from(prob_deletion_extend_artifact),
            indel_haplotype_window,
        )))
        .snv_read_evidence(RefCell::new(evidence::reads::SNVEvidence::new()))
        .mnv_read_evidence(RefCell::new(evidence::reads::MNVEvidence::new()))
        .indel_fragment_evidence(RefCell::new(evidence::fragments::IndelEvidence::new()))
        .none_read_evidence(RefCell::new(evidence::reads::NoneEvidence::new()))
    }
}

fn is_valid_record(record: &bam::Record) -> bool {
    !(record.is_secondary()
        || record.is_duplicate()
        || record.is_unmapped()
        || record.is_quality_check_failed())
}

impl Sample {
    /// Extract observations for the given variant.
    pub fn extract_observations(
        &mut self,
        start: u32,
        variant: &Variant,
        chrom: &[u8],
        chrom_seq: &[u8],
    ) -> Result<Pileup, Box<dyn Error>> {
        let centerpoint = variant.centerpoint(start);

        for vartype in &self.omit_repeat_regions {
            if variant.is_type(vartype) && is_repeat_variant(start, variant, chrom_seq) {
                // Do not return evidence, in order to mark variant in output as unclear.
                return Ok(Vec::new());
            }
        }

        self.record_buffer.fetch(
            chrom,
            start.saturating_sub(self.buffer_window),
            variant.end(start) + self.buffer_window,
        )?;

        let mut observations = Vec::new();

        match variant {
            //TODO: make &Variant::None add reads with position deleted if we want to check against indel alt alleles
            &Variant::SNV(_) | &Variant::MNV(_) | &Variant::None => {
                let mut candidate_records = Vec::new();
                // iterate over records
                for record in self.record_buffer.iter() {
                    if !is_valid_record(record) {
                        continue;
                    }
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
                    if !is_valid_record(record) {
                        continue;
                    }

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
                        if candidate.left.mapq() == 0 || right.mapq() == 0 {
                            // Ignore pairs with ambiguous alignments.
                            // The statistical model does not consider them anyway.
                            continue;
                        }
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
                            cigar.end_pos() as u32 + evidence::Clips::trailing(cigar).soft();

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
            }
        }

        // simulate strand bias
        // let mut count = 0;
        // for i in 0..observations.len() {
        //     if observations[i].prob_alt > observations[i].prob_ref && count < 6 {
        //         observations[i].reverse_strand = true;
        //         observations[i].forward_strand = false;
        //         count += 1;
        //     }
        // }

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
    ) -> Result<Option<Observation>, Box<dyn Error>> {
        let mut evidence: RefMut<dyn evidence::reads::AbstractReadEvidence> = match variant {
            &Variant::Deletion(_) | &Variant::Insertion(_) => self.indel_read_evidence.borrow_mut(),
            &Variant::SNV(_) => self.snv_read_evidence.borrow_mut(),
            &Variant::MNV(_) => self.mnv_read_evidence.borrow_mut(),
            &Variant::None => self.none_read_evidence.borrow_mut(),
        };

        if let Some((prob_ref, prob_alt)) =
            evidence.prob(record, cigar, start, variant, chrom_seq)?
        {
            let (prob_mapping, _) = evidence.prob_mapping_mismapping(record);

            // METHOD: This is an estimate of the allele likelihood at the true location in case
            // the read is mismapped. The value has to be approximately in the range of prob_alt
            // and prob_ref. Otherwise it could cause numerical problems, by dominating the
            // likelihood such that subtle differences in allele frequencies become numercically
            // invisible in the resulting likelihood.
            let prob_missed_allele = prob_ref.ln_add_exp(prob_alt) - LogProb(2.0_f64.ln());

            let prob_sample_alt = evidence.prob_sample_alt(
                record.seq().len() as u32,
                variant,
                &self.alignment_properties,
            );
            let strand = evidence.strand(record);
            dbg!(&strand);
            Ok(Some(
                ObservationBuilder::default()
                    .prob_mapping_mismapping(prob_mapping)
                    .prob_alt(prob_alt)
                    .prob_ref(prob_ref)
                    .prob_missed_allele(prob_missed_allele)
                    .prob_sample_alt(prob_sample_alt)
                    .prob_overlap(LogProb::ln_zero()) // no double overlap possible
                    .prob_any_strand(LogProb::from(Prob(0.5)))
                    .forward_strand(strand == Strand::Forward)
                    .reverse_strand(strand == Strand::Reverse)
                    .build()
                    .unwrap(),
            ))
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
    ) -> Result<Option<Observation>, Box<dyn Error>> {
        let prob_read = |record: &bam::Record,
                         cigar: &CigarStringView|
         -> Result<(LogProb, LogProb), Box<dyn Error>> {
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

        // METHOD: This is an estimate of the allele likelihood at the true location in case
        // the read is mismapped. The value has to be approximately in the range of prob_alt
        // and prob_ref. Otherwise it could cause numerical problems, by dominating the
        // likelihood such that subtle differences in allele frequencies become numercically
        // invisible in the resulting likelihood.
        let p_missed_left = p_ref_left.ln_add_exp(p_alt_left) - LogProb(2.0_f64.ln());
        let p_missed_right = p_ref_right.ln_add_exp(p_alt_right) - LogProb(2.0_f64.ln());

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

        let mut left_strand = self.indel_read_evidence.borrow().strand(left_record);
        let mut right_strand = self.indel_read_evidence.borrow().strand(right_record);
        // TODO find a better way to detect if there was no relevant overlap
        if p_alt_left == p_ref_left {
            left_strand = Strand::Unknown;
        }
        if p_alt_right == p_ref_right {
            right_strand = Strand::Unknown;
        }
        let mut forward_strand = left_strand == Strand::Forward || right_strand == Strand::Forward;
        let mut reverse_strand = left_strand == Strand::Reverse || right_strand == Strand::Reverse;
        if !forward_strand && !reverse_strand {
            // If there is no stranded evidence at all, consider observation to come from any
            // of the two strands.
            forward_strand = true;
            reverse_strand = true;
        }
        let prob_double_overlap = self.indel_fragment_evidence.borrow().prob_double_overlap(
            left_read_len,
            right_read_len,
            variant,
            &self.alignment_properties,
        );

        let prob_any_strand = LogProb::from(match self.protocol_strandedness {
            ProtocolStrandedness::Opposite => Prob(1.0 / 3.0),
            ProtocolStrandedness::Same => Prob(0.5),
        });

        let obs = ObservationBuilder::default()
            .prob_mapping_mismapping(prob_mismapping.ln_one_minus_exp())
            .prob_alt(p_alt_isize + p_alt_left + p_alt_right)
            .prob_ref(p_ref_isize + p_ref_left + p_ref_right)
            .prob_missed_allele(p_missed_left + p_missed_right)
            .prob_sample_alt(prob_sample_alt)
            .prob_overlap(prob_double_overlap)
            .prob_any_strand(prob_any_strand)
            .forward_strand(forward_strand)
            .reverse_strand(reverse_strand)
            .build()
            .unwrap();

        assert!(obs.prob_alt.is_valid());
        assert!(obs.prob_ref.is_valid());

        Ok(Some(obs))
    }
}

#[cfg(test)]
mod tests {
    extern crate env_logger;

    use super::*;
    use crate::constants;
    use crate::model;

    use crate::estimation::alignment_properties::{AlignmentProperties, InsertSize};
    use bio::io::fasta::{self, FastaRead};
    use bio::stats::{LogProb, PHREDProb, Prob};
    use itertools::Itertools;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::str;

    fn setup_sample(isize_mean: f64) -> Sample {
        SampleBuilder::default()
            .alignments(
                bam::IndexedReader::from_path(&"tests/indels.bam").unwrap(),
                AlignmentProperties::default(InsertSize {
                    mean: isize_mean,
                    sd: 20.0,
                }),
            )
            .error_probs(
                constants::PROB_ILLUMINA_INS,
                constants::PROB_ILLUMINA_DEL,
                Prob(0.0),
                Prob(0.0),
                100,
            )
            .max_depth(200)
            .build()
            .unwrap()
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

        for (record, true_alt_prob) in records.into_iter().zip(true_alt_probs.iter()) {
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

// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::f64;
use std::path::Path;
use std::str;
use std::hash::Hash;

use anyhow::Result;
use bio::stats::{LogProb, Prob};
use bio_types::{genome, genome::AbstractInterval};
use derive_builder::Builder;
use rand::distributions;
use rand::distributions::Distribution;
use rand::{rngs::StdRng, SeedableRng};
use rust_htslib::bam;

use crate::estimation::alignment_properties;
use crate::model::evidence;
use crate::model::evidence::{Observation, observation::Observable, observation};
use crate::model::VariantType;
use crate::variants::{Variant, self};

#[derive(new, Getters, Debug)]
pub struct RecordBuffer {
    inner: bam::RecordBuffer,
    #[getset(get = "pub")]
    single_read_window: u64,
    #[getset(get = "pub")]
    read_pair_window: u64,
}

impl RecordBuffer {
    pub fn window(&self, read_pair_mode: bool) -> u64 {
        if read_pair_mode { self.read_pair_window } else {self.single_read_window}
    }

    pub fn fetch(&mut self, interval: &genome::Interval, read_pair_mode: bool) -> Result<()> {
        self.inner.fetch(
            interval.contig().as_bytes(),
            interval.range().start.saturating_sub(self.window(read_pair_mode)),
            interval.range().end + self.window(read_pair_mode),
        )?;

        Ok(())
    }

    pub fn build_fetches(&self, read_pair_mode: bool) -> Fetches {
        Fetches { fetches: Vec::new(), window: self.window(read_pair_mode) }
    }

    pub fn iter(&self) -> impl Iterator<Item = &bam::Record> {
        self.inner.iter().filter(|record| is_valid_record(record))
    }
}

#[derive(Default, Derefable)]
pub struct Fetches {
    #[deref]
    fetches: Vec<genome::Interval>,
    window: u64
}

impl Fetches {
    pub fn push(&mut self, interval: &genome::Interval) {
        if let Some(last) = self.fetches.last_mut() {
            if last.contig() == interval.contig() {
                if interval.range().start.saturating_sub(self.window) <= last.range().end + self.window {
                    // merge the two intervals
                    last.range_mut().end = interval.range().end;
                    return;
                }
            }
        }

        self.fetches.push(interval.to_owned());
    }
}

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
) -> Result<alignment_properties::AlignmentProperties> {
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
    record_buffer: RecordBuffer,
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
        let read_pair_window = (alignment_properties.insert_size().mean
            + alignment_properties.insert_size().sd * 6.0) as u64;
        let single_read_window = alignment_properties.max_read_len as u64;
        self.alignment_properties(alignment_properties)
            .record_buffer(RecordBuffer::new(bam::RecordBuffer::new(bam, true), single_read_window, read_pair_window))
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
    pub fn extract_observations<'a, V, E, L>(
        &'a mut self,
        variant: &V,
    ) -> Result<Pileup>
    where 
        E: observation::Evidence<'a> + Eq + Hash,
        L: variants::Loci,
        V: Variant<'a, Loci=L, Evidence=E> + Observable<'a, E>,
    {
        variant.extract_observations(
            &mut self.record_buffer,
            &mut self.alignment_properties,
            self.max_depth,
        )
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
                AlignmentProperties::default(Some(InsertSize {
                    mean: isize_mean,
                    sd: 20.0,
                })),
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

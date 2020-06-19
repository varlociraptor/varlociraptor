// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::f64;
use std::hash::Hash;
use std::path::Path;
use std::rc::Rc;
use std::str;

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
use crate::model::evidence::{observation, observation::Observable, Observation};
use crate::model::VariantType;
use crate::variants::{self, Variant};

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
        if read_pair_mode {
            self.read_pair_window
        } else {
            self.single_read_window
        }
    }

    pub fn fetch(&mut self, interval: &genome::Interval, read_pair_mode: bool) -> Result<()> {
        dbg!("interval");
        dbg!(interval
            .range()
            .start
            .saturating_sub(self.window(read_pair_mode)));
        dbg!(interval.range().end + self.window(read_pair_mode));
        self.inner.fetch(
            interval.contig().as_bytes(),
            interval
                .range()
                .start
                .saturating_sub(self.window(read_pair_mode)),
            interval.range().end + self.window(read_pair_mode),
        )?;

        Ok(())
    }

    pub fn build_fetches(&self, read_pair_mode: bool) -> Fetches {
        Fetches {
            fetches: Vec::new(),
            window: self.window(read_pair_mode),
        }
    }

    pub fn iter<'a>(&'a self) -> impl Iterator<Item = Rc<bam::Record>> + 'a {
        self.inner
            .iter()
            .filter(|record| is_valid_record(record.as_ref()))
            .map(|record| Rc::clone(record))
    }
}

#[derive(Default, Derefable)]
pub struct Fetches {
    #[deref]
    fetches: Vec<genome::Interval>,
    window: u64,
}

impl Fetches {
    pub fn push(&mut self, interval: &genome::Interval) {
        if let Some(last) = self.fetches.last_mut() {
            if last.contig() == interval.contig() {
                if interval.range().start.saturating_sub(self.window)
                    <= last.range().end + self.window
                {
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
    Display,
    Debug,
    Clone,
    Copy,
    Serialize,
    Deserialize,
    EnumString,
    EnumIter,
    IntoStaticStr,
    EnumVariantNames,
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
    omit_insert_size: bool,
) -> Result<alignment_properties::AlignmentProperties> {
    let mut bam = bam::Reader::from_path(path)?;
    Ok(alignment_properties::AlignmentProperties::estimate(
        &mut bam,
        omit_insert_size,
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
            .record_buffer(RecordBuffer::new(
                bam::RecordBuffer::new(bam, true),
                single_read_window,
                read_pair_window,
            ))
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
    pub fn extract_observations<V, E, L>(&mut self, variant: &V) -> Result<Pileup>
    where
        E: observation::Evidence + Eq + Hash,
        L: variants::Loci,
        V: Variant<Loci = L, Evidence = E> + Observable<E>,
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
            .max_depth(200)
            .build()
            .unwrap()
    }

    fn ref_seq() -> Vec<u8> {
        let mut fa = fasta::Reader::from_file(&"tests/chr17.prefix.fa").unwrap();
        let mut chr17 = fasta::Record::new();
        fa.read(&mut chr17).unwrap();

        chr17.seq().to_owned()
    }
}

// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use super::evidence::observations::fragment_id_factory::FragmentIdFactory;
use super::evidence::observations::read_observation::major_alt_locus;
use super::evidence::realignment::Realignable;
use super::types::methylation::{meth_pos, meth_probs};
use crate::estimation::alignment_properties;
use crate::reference;
use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::{
    self, major_read_position, Observable, ReadObservation,
};
use crate::variants::types::Loci;
use crate::variants::{self, types::Variant};
use anyhow::Result;
use bio::stats::LogProb;
use bio_types::{genome, genome::AbstractInterval};
use by_address::ByAddress;
use derive_builder::Builder;
use rand::distributions;
use rand::distributions::Distribution;
use rand::{rngs::StdRng, SeedableRng};
use rust_htslib::bam::{self, Record};
use std::collections::{HashMap, HashSet};
use std::f64;
use std::path::Path;
use std::rc::Rc;
use std::str;

type MethylationPosToProbs = HashMap<usize, LogProb>;
type MethylationOfRead = HashMap<ByAddress<Rc<Record>>, Option<Rc<MethylationPosToProbs>>>;

#[derive(Getters, Debug)]
pub(crate) struct RecordBuffer {
    inner: bam::RecordBuffer,
    #[getset(get = "pub")]
    single_read_window: u64,
    #[getset(get = "pub")]
    read_pair_window: u64,
    #[getset(get = "pub")]
    methylation_probs: Option<MethylationOfRead>,
    #[getset(get = "pub")]
    failed_reads: Option<Vec<ByAddress<Rc<Record>>>>,
}

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct RecId {
    qname: String,
    pos: i64,
    flag: u16,
    tid: i32,
    // cigar: String,
}

impl RecId {}

impl RecordBuffer {
    pub(crate) fn new(
        inner: bam::RecordBuffer,
        single_read_window: u64,
        read_pair_window: u64,
        methylation_mm_ml_tag: bool,
    ) -> Self {
        RecordBuffer {
            inner,
            single_read_window,
            read_pair_window,
            methylation_probs: if methylation_mm_ml_tag {
                Some(HashMap::new())
            } else {
                None
            },
            failed_reads: if methylation_mm_ml_tag {
                Some(Vec::new())
            } else {
                None
            },
        }
    }

    pub(crate) fn window(&self, read_pair_mode: bool, left: bool) -> u64 {
        if read_pair_mode {
            self.read_pair_window
        } else if left {
            self.single_read_window
        } else {
            0
        }
    }

    pub(crate) fn get_methylation_probs(
        &self,
        rec: Rc<Record>,
    ) -> Option<Rc<HashMap<usize, LogProb>>> {
        self.methylation_probs
            .as_ref()
            .and_then(|meth_probs| meth_probs.get(&ByAddress(rec)).cloned().flatten())
    }

    pub(crate) fn fetch(
        &mut self,
        interval: &genome::Interval,
        read_pair_mode: bool,
    ) -> Result<()> {
        self.inner.fetch(
            interval.contig().as_bytes(),
            interval
                .range()
                .start
                .saturating_sub(self.window(read_pair_mode, true)),
            interval.range().end + self.window(read_pair_mode, false),
        )?;
        // If we are interested in methylation on PacBio or Nanopore data we need to compute the methylation probabilities
        if let Some(methylation_probs) = &mut self.methylation_probs {
            if let Some(failed_reads) = &mut self.failed_reads {
                // let mut first_it = true;
                let mut rec_debug = Vec::new();
                for rec in self.inner.iter() {
                    let rec_id = ByAddress(rec.clone());
                    rec_debug.push(rec.inner.core.pos);
                    // Compute methylation probs out of MM and ML tag and save in methylation_probs
                    if methylation_probs.get(&rec_id).is_none() && !failed_reads.contains(&rec_id) {
                        let pos_to_probs = meth_pos(rec).and_then(|meth_pos| {
                            meth_probs(rec).map(|meth_probs| {
                                let map: HashMap<usize, LogProb> =
                                    meth_pos.into_iter().zip(meth_probs.into_iter()).collect();
                                Rc::new(map) // Wrap the HashMap in Rc
                            })
                        });

                        if pos_to_probs.is_none() {
                            failed_reads.push(rec_id);
                        } else {
                            methylation_probs.insert(rec_id, pos_to_probs);
                        }
                    }
                }
                let buffer_ids: HashSet<_> = self
                    .inner
                    .iter()
                    .map(|rec| ByAddress(rec.clone()))
                    .collect();
                if let Some(methylation_probs_map) = &mut self.methylation_probs {
                    methylation_probs_map.retain(|key, _value| buffer_ids.contains(key));
                }
            }
        }

        Ok(())
    }

    pub(crate) fn build_fetches(&self, read_pair_mode: bool) -> Fetches {
        Fetches {
            fetches: Vec::new(),
            window: self.window(read_pair_mode, true),
        }
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = Rc<bam::Record>> + '_ {
        self.inner
            .iter()
            .filter(|record| is_valid_record(record.as_ref()))
            .map(Rc::clone)
    }
}

#[derive(Default, Derefable)]
pub(crate) struct Fetches {
    #[deref]
    fetches: Vec<genome::Interval>,
    window: u64,
}

impl Fetches {
    pub(crate) fn push(&mut self, interval: &genome::Interval) {
        if let Some(last) = self.fetches.last_mut() {
            if last.contig() == interval.contig()
                && interval.range().start.saturating_sub(self.window)
                    <= last.range().end + self.window
            {
                // merge the two intervals
                last.range_mut().end = interval.range().end;
                return;
            }
        }

        self.fetches.push(interval.to_owned());
    }
}

pub(crate) enum SubsampleCandidates {
    Necessary {
        rng: StdRng,
        prob: f64,
        prob_range: distributions::Uniform<f64>,
    },
    None,
}

impl SubsampleCandidates {
    pub(crate) fn new(max_depth: usize, depth: usize) -> Self {
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

    pub(crate) fn keep(&mut self) -> bool {
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

pub(crate) fn estimate_alignment_properties<P: AsRef<Path>>(
    paths: &[P],
    omit_insert_size: bool,
    reference_buffer: &mut reference::Buffer,
    num_records: Option<usize>,
) -> Result<alignment_properties::AlignmentProperties> {
    alignment_properties::AlignmentProperties::estimate(
        paths,
        omit_insert_size,
        reference_buffer,
        num_records,
    )
}

/// A sequenced sample, e.g., a tumor or a normal sample.
#[derive(Builder, Debug)]
#[builder(pattern = "owned")]
pub(crate) struct Sample {
    #[builder(private)]
    record_buffer: RecordBuffer,
    #[builder(private)]
    alignment_properties: alignment_properties::AlignmentProperties,
    #[builder(default = "200")]
    max_depth: usize,
    #[builder(default)]
    fragment_id_factory: FragmentIdFactory,
    report_fragment_ids: bool,
    adjust_prob_mapping: bool,
}

impl SampleBuilder {
    /// Register alignment information.
    ///
    /// # Arguments
    /// * `bam` - BAM file with the aligned and deduplicated sequence reads.
    pub(crate) fn alignments(
        self,
        bam: bam::IndexedReader,
        alignment_properties: alignment_properties::AlignmentProperties,
        min_refetch_distance: u64,
        methylation_mm_ml_tag: bool,
    ) -> Self {
        // METHOD: add maximum deletion len as this can make the footprint of the read on the reference
        // effectively larger. Additionally add some 10 bases further to account for uncertainty in the
        // estimated maximum deletion len.
        let single_read_window = alignment_properties.max_read_len as u64
            + alignment_properties
                .max_del_cigar_len
                .map_or(0, |l| l as u64)
            + 10;

        let read_pair_window = match alignment_properties.insert_size {
            Some(isize) => (isize.mean + isize.sd * 6.0) as u64,
            None => single_read_window,
        };
        let mut record_buffer = bam::RecordBuffer::new(bam, true);
        record_buffer.set_min_refetch_distance(min_refetch_distance);
        self.alignment_properties(alignment_properties)
            .record_buffer(RecordBuffer::new(
                record_buffer,
                single_read_window,
                read_pair_window,
                methylation_mm_ml_tag,
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
    pub(crate) fn extract_observations<V>(
        &mut self,
        variant: &V,
        alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Pileup>
    where
        V: Variant + Observable,
    {
        let mut observation_id_factory = if let Some(contig) = variant.loci().contig() {
            if self.report_fragment_ids {
                // METHOD: we only report read IDs for single contig variants.
                // Reason: we expect those to come in sorted, so that we can clear the
                // read ID registry at each new contig, saving lots of memory.
                // TODO: In the future, we might find a smarter way and thereby also include
                // multi-contig variants into the calculation.
                self.fragment_id_factory.register_contig(contig);
                Some(&mut self.fragment_id_factory)
            } else {
                None
            }
        } else {
            None
        };
        let observations = variant.extract_observations(
            &mut self.record_buffer,
            &mut self.alignment_properties,
            self.max_depth,
            alt_variants,
            &mut observation_id_factory,
        )?;
        // Process for each observation whether it is from the major read position or not.
        let major_pos = major_read_position(&observations);
        let major_alt_locus = major_alt_locus(&observations, &self.alignment_properties);
        let mut observations: Vec<_> = observations
            .iter()
            .map(|obs| obs.process(major_pos, &major_alt_locus, &self.alignment_properties))
            .collect();
        if self.adjust_prob_mapping {
            ReadObservation::adjust_prob_mapping(&mut observations, &self.alignment_properties);
        }
        Ok(Pileup::new(observations, Vec::new())) // TODO add depth observations!
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
pub enum Readtype {
    #[strum(serialize = "Nanopore")]
    Nanopore,
    #[strum(serialize = "Illumina")]
    Illumina,
    #[strum(serialize = "PacBio")]
    PacBio,
}

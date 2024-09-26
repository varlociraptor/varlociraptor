// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::f64;
use std::path::Path;
use std::rc::Rc;
use std::str;

use anyhow::Result;
use bio::io::om::bnx;
use bio::io::om::xmap;
use bio_types::{genome, genome::AbstractInterval};
use derive_builder::Builder;
use rand::distributions;
use rand::distributions::Distribution;
use rand::{rngs::StdRng, SeedableRng};
use rust_htslib::bam;

use crate::estimation::alignment_properties;
use crate::reference;
use crate::variants::evidence::observations::read_observation::{
    major_read_position, Observable, ReadObservation,
};
use crate::variants::types::Variant;

use super::evidence::observations::fragment_id_factory::FragmentIdFactory;
use super::evidence::observations::read_observation::major_alt_locus;
use super::evidence::realignment::Realignable;
use super::types::Loci;
use crate::variants::evidence::observations::pileup::Pileup;

#[derive(new, Getters, Debug)]
pub(crate) struct SequencingRecordBuffer {
    inner: bam::RecordBuffer,
    #[getset(get = "pub")]
    single_read_window: u64,
    #[getset(get = "pub")]
    read_pair_window: u64,
}

impl SequencingRecordBuffer {
    pub(crate) fn window(&self, read_pair_mode: bool, left: bool) -> u64 {
        if read_pair_mode {
            self.read_pair_window
        } else if left {
            self.single_read_window
        } else {
            0
        }
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

#[derive(new, Getters, Debug)]
pub(crate) struct OpticalMappingRecordBuffer {
    xmap: xmap::Container,
    bnx: bnx::Container,
}

impl OpticalMappingRecordBuffer {
    pub(crate) fn fetch(
        &self,
        interval: &genome::Interval,
    ) -> Result<impl Iterator<Item = &Rc<xmap::Record>>> {
        self.xmap.fetch(
            interval.contig().parse::<u32>()?,
            interval.range().start,
            interval.range().end,
        )
    }

    pub(crate) fn qry_record(&self, bnx_id: &u32) -> Result<&Rc<bnx::Record>> {
        self.bnx.record(bnx_id.clone())
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
    path: P,
    omit_insert_size: bool,
    reference_buffer: &mut reference::Buffer,
    num_records: Option<usize>,
) -> Result<alignment_properties::AlignmentProperties> {
    alignment_properties::AlignmentProperties::estimate(
        path,
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
    record_buffer: SequencingRecordBuffer,
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
            .record_buffer(SequencingRecordBuffer::new(
                record_buffer,
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

        let observations = variant.extract_sequencing_read_observations(
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

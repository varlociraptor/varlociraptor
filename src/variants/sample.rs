// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::f64;
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::str;

use anyhow::Result;
use bio_types::{genome, genome::AbstractInterval};
use derive_builder::Builder;
use rand::distributions;
use rand::distributions::Distribution;
use rand::{rngs::StdRng, SeedableRng};
use rust_htslib::bam;

use crate::estimation::alignment_properties;
use crate::reference;
use crate::variants::evidence::observations::depth_observation::{
    DepthObservable, DepthObservation,
};
use crate::variants::evidence::observations::observation::{AltLocus, ReadPosition};
use crate::variants::evidence::observations::read_observation::{
    self, major_read_position as read_major_read_position, ReadObservable, ReadObservation,
};
use crate::variants::types::{DepthVariant, ReadVariant};

use super::evidence::observations::fragment_id_factory::FragmentIdFactory;
use super::evidence::observations::read_observation::major_alt_locus as read_major_alt_locus;
use super::evidence::realignment::Realignable;
use super::types::Loci;
use crate::variants::evidence::observations::pileup::Pileup;

#[derive(new, Getters, Debug)]
pub(crate) struct RecordBuffer {
    inner: bam::RecordBuffer,
    #[getset(get = "pub")]
    single_read_window: u64,
    #[getset(get = "pub")]
    read_pair_window: u64,
}

impl RecordBuffer {
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
    bam_paths: &[P],
    cnv_path: Option<P>,
    omit_insert_size: bool,
    reference_buffer: &mut reference::Buffer,
    num_records: Option<usize>,
) -> Result<alignment_properties::AlignmentProperties> {
    alignment_properties::AlignmentProperties::estimate(
        bam_paths,
        cnv_path,
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
    #[builder(private)]
    bam_path: std::path::PathBuf,
    #[builder(private)]
    min_refetch_distance: u64,
}

impl SampleBuilder {
    /// Register alignment information.
    ///
    /// # Arguments
    /// * `bam` - BAM file with the aligned and deduplicated sequence reads.
    pub(crate) fn alignments(
        self,
        bam_path: &PathBuf,
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
            .record_buffer(RecordBuffer::new(
                record_buffer,
                single_read_window,
                read_pair_window,
            ))
            .bam_path(bam_path.into()) // Pfad speichern
            .min_refetch_distance(min_refetch_distance)
    }
}

fn is_valid_record(record: &bam::Record) -> bool {
    !(record.is_secondary()
        || record.is_duplicate()
        || record.is_unmapped()
        || record.is_quality_check_failed())
}

impl Sample {
    fn combine_pileup(
        read_obs: Option<Vec<ReadObservation<ReadPosition, AltLocus>>>,
        depth_obs: Option<Vec<DepthObservation>>,
    ) -> Result<Pileup> {
        Ok(Pileup::new(
            read_obs.unwrap_or_default(),
            depth_obs.unwrap_or_default(),
        ))
    }

    // This method is required when CNV probabilities are computed and subsequently recomputed for deletions,
    // as the probability computation alters the state of the record_buffer and it needs to be reset.

    pub(crate) fn reset_record_buffer(&mut self) -> Result<()> {
        let bam_reader = bam::IndexedReader::from_path(&self.bam_path)?;
        let mut new_inner = bam::RecordBuffer::new(bam_reader, true);
        new_inner.set_min_refetch_distance(self.min_refetch_distance);

        self.record_buffer = RecordBuffer::new(
            new_inner,
            self.record_buffer.single_read_window,
            self.record_buffer.read_pair_window,
        );
        Ok(())
    }

    pub(crate) fn extract_read_and_depth_observations<V>(
        &mut self,
        variant: &V,
        alt_variants: &[Box<dyn Realignable>],
        max_number_cn: usize,
    ) -> Result<Pileup>
    where
        V: ReadVariant + ReadObservable + DepthVariant + DepthObservable,
    {
        let read_obs = self.compute_read_observations(variant, alt_variants)?;
        self.reset_record_buffer()?;
        let depth_obs = self.compute_depth_observations(variant, alt_variants, max_number_cn)?;
        Self::combine_pileup(Some(read_obs), Some(depth_obs))
    }

    pub(crate) fn extract_read_observations<V>(
        &mut self,
        variant: &V,
        alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Pileup>
    where
        V: ReadVariant + ReadObservable,
    {
        let read_obs = self.compute_read_observations(variant, alt_variants)?;
        Self::combine_pileup(Some(read_obs), None)
    }

    pub(crate) fn extract_depth_observations<V>(
        &mut self,
        variant: &V,
        alt_variants: &[Box<dyn Realignable>],
        max_number_cn: usize,
    ) -> Result<Pileup>
    where
        V: DepthVariant + DepthObservable,
    {
        let depth_obs = self.compute_depth_observations(variant, alt_variants, max_number_cn)?;
        Self::combine_pileup(None, Some(depth_obs))
    }

    /// Extract observations for the given variant.
    pub(crate) fn compute_read_observations<V>(
        &mut self,
        variant: &V,
        alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Vec<ReadObservation<ReadPosition, AltLocus>>>
    where
        V: ReadVariant + ReadObservable,
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
        let major_pos = read_major_read_position(&observations);
        let major_alt_locus = read_major_alt_locus(&observations, &self.alignment_properties);
        let mut observations: Vec<_> = observations
            .iter()
            .map(|obs| obs.process(major_pos, &major_alt_locus, &self.alignment_properties))
            .collect();
        if self.adjust_prob_mapping {
            ReadObservation::adjust_prob_mapping(&mut observations, &self.alignment_properties);
        }
        Ok(observations) // TODO add depth observations!
    }

    pub(crate) fn compute_depth_observations<V>(
        &mut self,
        variant: &V,
        alt_variants: &[Box<dyn Realignable>],
        max_number_cn: usize,
    ) -> Result<Vec<DepthObservation>>
    where
        V: DepthVariant + DepthObservable,
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
            max_number_cn,
        )?;
        // Process for each observation whether it is from the major read position or not.
        // let major_pos = depth_major_read_position(&observations);
        // let major_alt_locus = depth_major_alt_locus(&observations, &self.alignment_properties);
        // let mut observations: Vec<_> = observations
        //     .iter()
        //     .map(|obs| obs.process(major_pos, &major_alt_locus, &self.alignment_properties))
        //     .collect();
        // if self.adjust_prob_mapping {
        //     DepthObservation::adjust_prob_mapping(&mut observations, &self.alignment_properties);
        // }
        Ok(observations) // TODO add depth observations!
    }
}

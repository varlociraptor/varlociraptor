use std::collections::HashMap;
use std::path::PathBuf;
use std::rc::Rc;
use std::str;
use std::sync::RwLock;

use anyhow::{Context, Result};
use bio::stats::{bayesian, LogProb};
use bio_types::genome;
use bio_types::genome::AbstractLocus;
use derive_builder::Builder;
use itertools::Itertools;
use progress_logger::ProgressLogger;
use rust_htslib::bcf::{self, Read};

use crate::calling::variants::preprocessing::{
    read_observations, remove_observation_header_entries, OBSERVATION_FORMAT_VERSION,
};
use crate::calling::variants::SampleInfo;
use crate::calling::variants::{
    chrom, event_tag_name, Call, CallBuilder, SampleInfoBuilder, VariantBuilder,
};
use crate::errors;
use crate::grammar;
use crate::utils;
use crate::variants::evidence::observations::pileup::Pileup;

use crate::variants::model::modes::generic::LikelihoodOperands;
use crate::variants::model::modes::generic::{
    self, GenericLikelihood, GenericModelBuilder, GenericPosterior,
};
use crate::variants::model::Contamination;
use crate::variants::model::{self};
use crate::variants::model::{bias::Artifacts, AlleleFreq};
use crate::variants::types::breakends::BreakendIndex;

use super::preprocessing::Observations;

pub(crate) type AlleleFreqCombination = LikelihoodOperands;

pub(crate) type Model<Pr> =
    bayesian::Model<GenericLikelihood, Pr, GenericPosterior, generic::Cache>;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct Caller<Pr>
where
    Pr: bayesian::model::Prior,
{
    samplenames: grammar::SampleInfo<String>,
    observations: grammar::SampleInfo<Option<PathBuf>>,
    omit_strand_bias: bool,
    omit_read_orientation_bias: bool,
    omit_read_position_bias: bool,
    omit_softclip_bias: bool,
    omit_homopolymer_artifact_detection: bool,
    scenario: grammar::Scenario,
    outbcf: Option<PathBuf>,
    contaminations: grammar::SampleInfo<Option<Contamination>>,
    resolutions: grammar::SampleInfo<grammar::Resolution>,
    prior: Pr,
    breakend_index: BreakendIndex,
    #[builder(default)]
    breakend_results: RwLock<HashMap<Vec<u8>, BreakendResult>>,
    log_each_record: bool,
}

impl<Pr> Caller<Pr>
where
    Pr: bayesian::model::Prior<Event = AlleleFreqCombination>
        + model::prior::UpdatablePrior
        + model::prior::CheckablePrior
        + Clone
        + Default,
{
    pub(crate) fn n_samples(&self) -> usize {
        self.samplenames.len()
    }

    pub(crate) fn header(&self) -> Result<bcf::Header> {
        let mut header = bcf::Header::from_template(
            bcf::Reader::from_path(self.observations.first_not_none().as_ref().unwrap())?.header(),
        );

        remove_observation_header_entries(&mut header);

        // register samples
        for sample_name in self.samplenames.iter() {
            header.push_sample(sample_name.as_bytes());
        }

        // register events
        for event in self.scenario.events().keys() {
            header.push_record(
                format!(
                    "##INFO=<ID={},Number=A,Type=Float,\
                     Description=\"Posterior probability for event {} (PHRED)\">",
                    event_tag_name(event),
                    event
                )
                .as_bytes(),
            );
        }
        header.push_record(
            b"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,\
             Description=\"Posterior probability for any artifact, indicated by strand, read position, \
             read orientation, softclip bias, or divindel bias (PHRED). See the bias specific records below for \
             an explanation for each type of bias.\">",
        );
        header.push_record(
            b"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,\
             Description=\"Posterior probability for not having a variant (PHRED)\">",
        );

        // register sample specific tags
        header.push_record(
            b"##FORMAT=<ID=DP,Number=1,Type=Integer,\
              Description=\"Expected sequencing depth, while considering mapping uncertainty\">",
        );
        header.push_record(
            b"##FORMAT=<ID=AF,Number=A,Type=Float,\
              Description=\"Maximum a posteriori probability estimate of allele frequency\">",
        );
        header.push_record(
            b"##FORMAT=<ID=SOBS,Number=A,Type=String,\
              Description=\"Summary of simplified observations. Each entry is encoded as CB, with C being a count, \
              B being the posterior odds for the alt allele. \
              Posterior odds for alt allele of each fragment are given as extended Kass Raftery \
              scores: N=none, E=equal, B=barely, P=positive, S=strong, V=very strong (lower case if \
              probability for correct mapping of fragment is <95%). Note that we extend Kass Raftery scores with \
              a term for equality between the evidence of the two alleles (E=equal).\">",
        );
        header.push_record(
            b"##FORMAT=<ID=OBS,Number=A,Type=String,\
              Description=\"Summary of observations. Each entry is encoded as CBTASOPXI, with C being a count, \
              B being the posterior odds for the alt allele (see below), T being the type of alignment, encoded \
              as s=single end and p=paired end, A denoting whether the observations also map to an alternative locus \
              (# = most found alternative locus, * = other locus, . = no locus), \
              S being the strand that supports the observation (+, -, or * for both), \
              O being the read orientation (> = F1R2, < = F2R1, * = unknown, ! = non standard, e.g. R1F2), \
              P being the read position (^ = most found read position, * = any other position or position is irrelevant), \
              X denoting whether the respective alignments entail a softclip ($ = softclip, . = no soft clip), and \
              I denoting indel operations in the respective alignments against the alt allele \
              (* = some indel, . = no indel or information irrelevant for variant type). \
              Posterior odds for alt allele of each fragment are given as extended Kass Raftery \
              scores: N=none, E=equal, B=barely, P=positive, S=strong, V=very strong (lower case if \
              probability for correct mapping of fragment does not correspond to the maximum reported value by the mapper \
              (for bwa, this is usually 60 in PHRED scale)). Note that we extend Kass Raftery scores with \
              a term for equality between the evidence of the two alleles (E=equal).\">",
        );
        header.push_record(
            b"##FORMAT=<ID=SB,Number=A,Type=String,\
              Description=\"Strand bias estimate: + indicates that ALT allele is associated with \
              forward strand, - indicates that ALT allele is associated with reverse strand, \
              . indicates no strand bias. Strand bias is indicative for systematic sequencing \
              errors. Probability for strand bias is captured by the ARTIFACT event (PROB_ARTIFACT).\">",
        );
        header.push_record(
            b"##FORMAT=<ID=ROB,Number=A,Type=String,\
              Description=\"Read orientation bias estimate: > indicates that ALT allele is associated with \
              F1R2 orientation, < indicates that ALT allele is associated with F2R1 orientation, \
              . indicates no read orientation bias. Read orientation bias is indicative of Guanin \
              oxidation artifacts. Probability for read orientation bias is captured by the ARTIFACT \
              event (PROB_ARTIFACT).\">",
        );
        header.push_record(
            b"##FORMAT=<ID=RPB,Number=A,Type=String,\
              Description=\"Read position bias estimate: ^ indicates that ALT allele is associated with \
              the most found read position, . indicates that there is no read position bias.
              Read position bias is indicative of systematic sequencing errors, e.g. in a specific cycle. \
              Probability for read orientation bias is captured by the ARTIFACT \
              event (PROB_ARTIFACT).\">",
        );
        header.push_record(
            b"##FORMAT=<ID=SCB,Number=A,Type=String,\
              Description=\"Softclip bias estimate: $ indicates that ALT allele is associated with \
              with softclips in the same alignment, . indicates that there is no softclip bias.
              Softclip bias is indicative of systematic alignment errors, cause by a part of the read \
              that does not properly align to the reference (and is thus soft clipped). Note that \
              softclips can also be caused by structural variants. However, structural variants on the \
              same haplotype as e.g. an SNV should not cause a softclip bias, because there will usually \
              still be reads that do not reach the SV, thereby providing evidence against a softclip \
              bias. Probability for softclip bias is captured by the ARTIFACT \
              event (PROB_ARTIFACT).\">",
        );
        header.push_record(
            b"##FORMAT=<ID=HE,Number=A,Type=String,\
              Description=\"Homopolymer error estimate: * indicates that ALT allele is associated with \
              with homopolymer indel operations of varying length, . indicates that there is no homopolymer error.
              Homopolymer error is indicative of systematic PCR amplification errors. \
              Probability for such homopolymer artifacts is captured by the ARTIFACT \
              event (PROB_ARTIFACT).\">",
        );
        header.push_record(
            b"##FORMAT=<ID=ALB,Number=A,Type=String,\
              Description=\"Alt locus bias estimate: * indicates that ALT allele is systematically associated \
              with either MAPQs smaller than the maximum MAPQ or a major alternative alignment (XA tag) \
              reported by the used read mapper. \
              This would be indicative of ALT reads actually coming from another locus (e.g. some repeat, \
              a homology, a distant variant allele, or a CNV). \
              Probability for mapping quality bias is captured by the ARTIFACT \
              event (PROB_ARTIFACT).\">",
        );
        header.push_record(
            b"##FORMAT=<ID=AFD,Number=A,Type=String,\
              Description=\"Sampled posterior probability densities of allele frequencies in PHRED scale \
              (the smaller the higher, with 0 being equal to an unscaled probability of 1). \
              In the discrete case (no somatic mutation rate or continuous universe in the scenario), \
              these can be seen as posterior probabilities. Note that densities can be greater than one.\">",
        );

        Ok(header)
    }

    pub(crate) fn writer(&self) -> Result<bcf::Writer> {
        let header = self.header();

        Ok(if let Some(ref path) = self.outbcf {
            bcf::Writer::from_path(path, header.as_ref().unwrap(), false, bcf::Format::Bcf)
                .context(format!("Unable to write BCF to {}.", path.display()))?
        } else {
            bcf::Writer::from_stdout(header.as_ref().unwrap(), false, bcf::Format::Bcf)
                .context("Unable to write BCF to STDOUT.")?
        })
    }

    fn model(&self) -> Model<Pr> {
        GenericModelBuilder::default()
            // TODO allow to define prior in the grammar
            .prior(self.prior.clone())
            .contaminations(self.contaminations.clone())
            .resolutions(self.resolutions.clone())
            .build()
            .unwrap()
    }

    fn observations(&self) -> Result<grammar::SampleInfo<Option<bcf::Reader>>> {
        let mut observations = grammar::SampleInfo::default();
        for path in self.observations.iter() {
            if let Some(path) = path {
                let mut reader = bcf::Reader::from_path(path)?;
                reader.set_threads(1)?;
                observations.push(Some(reader));
            } else {
                observations.push(None);
            }
        }
        Ok(observations)
    }

    pub(crate) fn call(&self) -> Result<()> {
        let mut observations = self.observations()?;
        let mut bcf_writer = self.writer()?;
        bcf_writer.set_threads(1)?;

        // Check observation format.
        for obs_reader in observations.iter_not_none() {
            let mut valid = false;
            for record in obs_reader.header().header_records() {
                if let bcf::HeaderRecord::Generic { key, value } = record {
                    if key == "varlociraptor_observation_format_version"
                        && value == OBSERVATION_FORMAT_VERSION
                    {
                        valid = true;
                    }
                }
            }
            if !valid {
                return Err(errors::Error::InvalidObservationFormat.into());
            }
        }

        // data structures
        // For SNVs and MNVs we need a special model as here read orientation bias and read position bias needs to be considered.
        let mut models = HashMap::new();
        let mut events = HashMap::new();
        let mut last_rids = HashMap::new();

        // process calls
        let mut i = 0;
        let mut progress_logger = ProgressLogger::builder()
            .with_items_name("records")
            .with_frequency(std::time::Duration::from_secs(20))
            .start();
        loop {
            let mut records =
                observations.map(|reader| reader.as_ref().map(|reader| reader.empty_record()));
            let mut eof = Vec::new();
            for item in observations.iter_mut().zip(records.iter_mut()) {
                if let (Some(reader), Some(record)) = item {
                    eof.push(match reader.read(record) {
                        None => true,
                        Some(res) => {
                            res?;
                            false
                        }
                    });
                }
            }

            if eof.iter().all(|v| *v) {
                progress_logger.stop();
                return Ok(());
            } else if !eof.iter().all(|v| !v) {
                // only some are EOF, this is an error
                return Err(errors::Error::InconsistentObservations.into());
            }

            // ensure that all observation BCFs contain exactly the same calls
            let first_record = records.first_not_none()?;
            let current_rid = first_record.rid();
            let current_pos = first_record.pos();
            let current_alleles = first_record.alleles();
            for record in records[1..].iter().flatten() {
                if record.rid() != current_rid
                    || record.pos() != current_pos
                    || record.alleles() != current_alleles
                {
                    return Err(errors::Error::InconsistentObservations.into());
                }
            }

            if self.log_each_record {
                info!(
                    "Processing record at {}:{}",
                    first_record.contig(),
                    first_record.pos() + 1
                );
            }

            // obtain variant type
            let variant_type =
                utils::collect_variants(records.first_not_none_mut()?, false, None)?[0].to_type();

            let mut work_item = self.preprocess_record(&mut records, i, &observations)?;

            // process work item
            let contig = str::from_utf8(work_item.call.chrom()).unwrap();
            let _model;
            let _events;
            let _last_rid;

            let model_mode = (
                work_item.check_read_orientation_bias,
                work_item.check_read_position_bias,
                work_item.check_softclip_bias,
                work_item.check_homopolymer_artifact_detection,
            );
            _model = models.entry(model_mode).or_insert_with(|| self.model());
            _events = events.entry(model_mode).or_insert_with(Vec::new);
            {
                let entry = last_rids.entry(model_mode).or_insert(None);
                _last_rid = *entry;
                *entry = Some(work_item.rid);
            }

            self.configure_model(
                work_item.rid,
                _last_rid,
                _model,
                _events,
                contig,
                variant_type,
                work_item.check_read_orientation_bias,
                work_item.check_strand_bias,
                work_item.check_read_position_bias,
                work_item.check_softclip_bias,
                work_item.check_homopolymer_artifact_detection,
            )?;

            self.call_record(&mut work_item, _model, &_events);

            work_item.call.write_final_record(&mut bcf_writer)?;
            progress_logger.update(1u64);

            i += 1;
        }
    }

    fn preprocess_record(
        &self,
        records: &mut grammar::SampleInfo<Option<bcf::Record>>,
        index: usize,
        observations: &grammar::SampleInfo<Option<bcf::Reader>>,
    ) -> Result<WorkItem> {
        let (call, snv, bnd_event, rid, is_snv_or_mnv) = {
            let first_record = records.first_not_none_mut()?;
            let start = first_record.pos() as u64;
            let chrom = chrom(observations.first_not_none()?, first_record);

            let _locus = genome::Locus::new(str::from_utf8(chrom).unwrap().to_owned(), start);

            let call = CallBuilder::default()
                .chrom(chrom.to_owned())
                .pos(start)
                .id({
                    let id = first_record.id();
                    if id == b"." {
                        None
                    } else {
                        Some(id)
                    }
                })
                .record(first_record)?
                .build()
                .unwrap();

            // store information about SNV for special handling in posterior computation (variant selection operations)
            let snv;
            let is_snv_or_mnv;
            {
                let alleles = first_record.alleles();
                if alleles[0].len() == 1 && alleles[1].len() == 1 {
                    is_snv_or_mnv = true;
                    snv = Some(model::modes::generic::Snv::new(
                        alleles[0][0],
                        alleles[1][0],
                    ));
                } else if alleles[0].len() == alleles[1].len() {
                    is_snv_or_mnv = true;
                    snv = None;
                } else {
                    is_snv_or_mnv = false;
                    snv = None;
                }
            };

            let bnd_event = if utils::is_bnd(first_record)? {
                Some(utils::info_tag_event(first_record)?.unwrap())
            } else {
                None
            };

            let rid = first_record
                .rid()
                .ok_or_else(|| errors::Error::RecordMissingChrom { i: index + 1 })?;

            (call, snv, bnd_event, rid, is_snv_or_mnv)
        };

        let mut variant_builder = VariantBuilder::default();
        variant_builder.record(records.first_not_none_mut()?)?;

        let mut work_item = WorkItem {
            rid,
            call,
            pileups: None,
            snv,
            bnd_event,
            variant_builder,
            index,
            check_read_orientation_bias: is_snv_or_mnv && !self.omit_read_orientation_bias,
            check_strand_bias: !self.omit_strand_bias,
            check_read_position_bias: is_snv_or_mnv && !self.omit_read_position_bias,
            check_softclip_bias: is_snv_or_mnv && !self.omit_softclip_bias,
            check_homopolymer_artifact_detection: false,
        };

        if let Some(ref event) = work_item.bnd_event {
            if self.breakend_results.read().unwrap().contains_key(event) {
                // METHOD: Another breakend in the same event was already processed, hence, we will just copy the
                // results (no pileup needed).
                return Ok(work_item);
            }
        }

        // obtain pileups
        let mut pileups = Vec::new();
        for record in records.iter_mut() {
            let pileup = if let Some(record) = record {
                let Observations {
                    mut pileup,
                    is_homopolymer_indel,
                } = read_observations(record)?;
                if is_homopolymer_indel && !self.omit_homopolymer_artifact_detection {
                    // METHOD: check for homopolymer artifacts if at least one pileup contains the corresponding information.
                    work_item.check_homopolymer_artifact_detection |= true;
                }
                if is_snv_or_mnv {
                    // METHOD: remove non-standard alignments. They might come from near
                    // SVs and can induce artifactual SNVs or MNVs. By removing them,
                    // we just conservatively reduce the coverage to those which are
                    // clearly not influenced by a close SV.
                    pileup.remove_nonstandard_alignments(self.omit_read_orientation_bias);
                }

                pileup
            } else {
                Pileup::default()
            };
            pileups.push(pileup);
        }

        work_item.pileups = Some(pileups);

        Ok(work_item)
    }

    fn configure_model(
        &self,
        current_rid: u32,
        rid: Option<u32>,
        model: &mut Model<Pr>,
        events: &mut Vec<model::Event>,
        contig: &str,
        variant_type: model::VariantType,
        consider_read_orientation_bias: bool,
        consider_strand_bias: bool,
        consider_read_position_bias: bool,
        consider_softclip_bias: bool,
        consider_homopolymer_error: bool,
    ) -> Result<()> {
        if !rid.map_or(false, |rid: u32| current_rid == rid) || events.is_empty() {
            // rid is not the same as before or the model mode has changed to something new, obtain event universe
            // clear old events
            events.clear();

            // register absent event
            events.push(model::Event {
                name: "absent".to_owned(),
                vafs: grammar::VAFTree::absent(self.n_samples()),
                biases: vec![Artifacts::none()],
            });

            // add events from scenario
            for (event_name, vaftree) in self.scenario.vaftrees(contig)? {
                events.push(model::Event {
                    name: event_name.clone(),
                    vafs: vaftree.clone(),
                    biases: vec![Artifacts::none()],
                });

                let biases: Vec<_> = Artifacts::all_artifact_combinations(
                    consider_read_orientation_bias,
                    consider_strand_bias,
                    consider_read_position_bias,
                    consider_softclip_bias,
                    consider_homopolymer_error,
                )
                .collect();

                if !biases.is_empty() {
                    // Corresponding biased event.
                    events.push(model::Event {
                        name: event_name.clone(),
                        vafs: vaftree.clone(),
                        biases,
                    });
                }
            }

            // update prior to the VAF universe of the current chromosome
            let mut vaf_universes = self.scenario.sample_info();
            let mut ploidies = self.scenario.sample_info();
            for (sample_name, sample) in self.scenario.samples().iter() {
                let universe = sample.contig_universe(&contig, self.scenario.species())?;
                vaf_universes = vaf_universes.push(sample_name, universe.to_owned());

                let ploidy = sample.contig_ploidy(&contig, self.scenario.species())?;
                ploidies = ploidies.push(sample_name, ploidy);
            }

            model
                .prior_mut()
                .set_universe_and_ploidies(vaf_universes.build(), ploidies.build());
            model.prior().check()?;
        }

        model.prior_mut().set_variant_type(variant_type);

        Ok(())
    }

    fn call_record(
        &self,
        work_item: &mut WorkItem,
        model: &Model<Pr>,
        event_universe: &[model::Event],
    ) {
        if let Some(ref bnd_event) = work_item.bnd_event {
            if let Some(result) = self.breakend_results.read().unwrap().get(bnd_event) {
                // Take sample info and event probs from previous breakend.
                work_item
                    .variant_builder
                    .event_probs(Some(result.event_probs.clone()));
                work_item
                    .variant_builder
                    .sample_info(result.sample_info.clone());

                let variant = work_item.variant_builder.build().unwrap();
                work_item.call.variant = Some(variant);

                return;
            }
        }

        if work_item.pileups.is_some() {
            let data = model::modes::generic::Data::new(
                work_item.pileups.take().unwrap(),
                work_item.snv.clone(),
            );

            let mut event_universe: Vec<_> = event_universe.to_vec();
            for event in &mut event_universe {
                // METHOD: learn parameters for each bias (if necessary).
                // By this, we can avoid marginalization of them, which is
                // unnecessarily expensive.
                for bias in &mut event.biases {
                    bias.learn_parameters(data.pileups());
                }
            }

            // Compute probabilities for given events.
            let m = model.compute(event_universe.iter().cloned(), &data);

            // add calling results
            let mut event_probs: HashMap<String, LogProb> = event_universe
                .iter()
                .filter_map(|event| {
                    if event.is_artifact() {
                        None
                    } else {
                        let p = m.posterior(event).unwrap();
                        Some((event.name.clone(), p))
                    }
                })
                .collect();

            // generate artifact event
            let prob_artifact = LogProb::ln_sum_exp(
                &event_universe
                    .iter()
                    .filter_map(|event| {
                        if event.is_artifact() {
                            let p = m.posterior(event).unwrap();
                            Some(p)
                        } else {
                            None
                        }
                    })
                    .collect_vec(),
            );

            event_probs.insert("artifact".to_owned(), prob_artifact);

            let is_artifact = event_probs
                .iter()
                .all(|(event, prob)| event == "artifact" || *prob < prob_artifact);

            work_item.variant_builder.event_probs(Some(event_probs));

            // add sample specific information
            work_item
                .variant_builder
                .sample_info(self.sample_infos(&m, is_artifact, data));
        } else {
            unreachable!();
        }

        let variant = work_item.variant_builder.build().unwrap();

        if let Some(ref event) = work_item.bnd_event {
            if self.breakend_index.last_record_index(event).unwrap() == work_item.index {
                // METHOD: last index, hence clear result
                self.breakend_results.write().unwrap().remove(event);
            } else {
                // METHOD: store breakend group result for next breakend of this group
                self.breakend_results.write().unwrap().insert(
                    event.to_owned(),
                    BreakendResult {
                        event_probs: variant.event_probs().as_ref().unwrap().clone(),
                        sample_info: variant.sample_info().clone(),
                    },
                );
            }
        }

        work_item.call.variant = Some(variant);
    }

    fn sample_infos(
        &self,
        model_instance: &bayesian::model::ModelInstance<AlleleFreqCombination, model::Event>,
        is_artifact: bool,
        data: model::modes::generic::Data,
    ) -> Vec<Option<SampleInfo>> {
        for (map_estimates, _) in model_instance.event_posteriors() {
            if map_estimates
                .iter()
                .any(|estimate| estimate.artifacts.is_artifact())
                && !is_artifact
            {
                // METHOD: skip MAP that is an artifact if the overall artifact event is not the strongest one.
                // This ensures consistency between the events and the per sample MAPs.
                continue;
            }
            // This is the MAP estimate we want to report, build sample information with it.
            return data
                .into_pileups()
                .into_iter()
                .zip(map_estimates.iter())
                .enumerate()
                .map(|(sample, (pileup, estimate))| {
                    let mut sample_builder = SampleInfoBuilder::default();
                    sample_builder.pileup(Rc::new(pileup));
                    match estimate {
                        model::likelihood::Event {
                            artifacts: biases, ..
                        } if biases.is_artifact() => {
                            sample_builder
                                .allelefreq_estimate(AlleleFreq(0.0))
                                .artifacts(biases.clone());
                        }
                        model::likelihood::Event { allele_freq, .. } => {
                            sample_builder
                                .allelefreq_estimate(*allele_freq)
                                .artifacts(Artifacts::none());
                        }
                    };

                    // METHOD: Collect VAF dist for sample by looking at all alternative VAFs with the same VAFs for the other samples as the MAP.
                    sample_builder.vaf_dist(if !estimate.is_artifact() {
                        Some(
                            model_instance
                                .event_posteriors()
                                .filter_map(|(estimate, prob)| {
                                    let event = estimate.events().get(sample).unwrap();
                                    let others_equal_map = || {
                                        map_estimates.events().iter().all(
                                            |(other_sample, map_event)| {
                                                //dbg!((other_sample, sample, estimate.events().get(other_sample).unwrap(), map_event));
                                                if other_sample != sample {
                                                    // check if other event is the same as the map event
                                                    let other_event = estimate
                                                        .events()
                                                        .get(other_sample)
                                                        .unwrap();
                                                    other_event == map_event
                                                } else {
                                                    // don't do that for our current sample
                                                    true
                                                }
                                            },
                                        )
                                    };
                                    if !event.is_artifact() && others_equal_map() {
                                        Some((event.allele_freq, prob))
                                    } else {
                                        None
                                    }
                                })
                                .collect(),
                        )
                    } else {
                        None
                    });

                    Some(sample_builder.build().unwrap())
                })
                .collect_vec();
        }

        // no observations
        vec![None; data.into_pileups().len()]
    }
}

#[derive(Default)]
pub(crate) struct BreakendResult {
    event_probs: HashMap<String, LogProb>,
    sample_info: Vec<Option<SampleInfo>>,
}

struct WorkItem {
    rid: u32,
    call: Call,
    variant_builder: VariantBuilder,
    pileups: Option<Vec<Pileup>>,
    snv: Option<model::modes::generic::Snv>,
    bnd_event: Option<Vec<u8>>,
    index: usize,
    check_read_orientation_bias: bool,
    check_strand_bias: bool,
    check_read_position_bias: bool,
    check_softclip_bias: bool,
    check_homopolymer_artifact_detection: bool,
}

// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

pub(crate) mod calling;
pub mod preprocessing;

use std::collections::HashMap;
use std::collections::HashSet;
use std::convert::TryFrom;
use std::ops::RangeInclusive;
use std::rc::Rc;
use std::str;

use anyhow::Result;
use bio::stats::{LogProb, PHREDProb};
use bio_types::sequence::SequenceReadPairOrientation;
use derive_builder::Builder;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use rust_htslib::bcf::{self, record::Numeric, Read};
use vec_map::VecMap;

use crate::calling::variants::preprocessing::write_observations;
use crate::utils;
use crate::utils::aux_info::AuxInfo;
use crate::utils::bayes_factor_to_letter;
use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::expected_depth;
use crate::variants::evidence::observations::read_observation::AltLocus;
use crate::variants::evidence::observations::read_observation::{ReadPosition, Strand};
use crate::variants::model;
use crate::variants::model::bias::AltLocusBias;
use crate::variants::model::HaplotypeIdentifier;
use crate::variants::model::VariantPrecision;
use crate::variants::model::{
    bias::Artifacts, bias::HomopolymerError, bias::ReadOrientationBias, bias::ReadPositionBias,
    bias::SoftclipBias, bias::StrandBias, AlleleFreq,
};

lazy_static! {
    static ref OMIT_AUX_INFO: HashSet<Vec<u8>> = HashSet::from([
        b"MATEID".to_vec(),
        b"EVENT".to_vec(),
        b"SVLEN".to_vec(),
        b"SVTYPE".to_vec(),
        b"END".to_vec(),
    ]);
}

#[derive(EnumString, Hash, Eq, PartialEq, Debug, Clone, Copy, IntoStaticStr, Display)]
#[strum(serialize_all = "kebab_case")]
pub(crate) enum Hint {
    AdjustedSingletonEvidence,
    FilteredNonStandardAlignments,
    MissingData,
}

#[derive(Default, Clone, Debug, Builder, Getters, CopyGetters)]
pub(crate) struct Call {
    #[getset(get = "pub(crate)")]
    chrom: Vec<u8>,
    #[getset(get = "pub(crate)")]
    pos: u64,
    #[builder(default = "None")]
    #[getset(get = "pub(crate)")]
    id: Option<Vec<u8>>,
    #[builder(default = "None")]
    #[getset(get = "pub(crate)")]
    mateid: Option<Vec<u8>>,
    #[builder(default = "None")]
    #[getset(get_copy = "pub(crate)")]
    heterozygosity: Option<LogProb>,
    #[builder(default = "None")]
    #[getset(get_copy = "pub(crate)")]
    somatic_effective_mutation_rate: Option<LogProb>,
    #[builder(default)]
    #[getset(get = "pub(crate)")]
    aux_info: AuxInfo,
    #[builder(default)]
    hints: HashSet<Hint>,
    //aux_fields: HashSet<Vec<u8>>,
    #[builder(default)]
    #[getset(get = "pub(crate)")]
    variant: Option<Variant>,
}

impl CallBuilder {
    pub(crate) fn record(&mut self, record: &mut bcf::Record) -> Result<&mut Self> {
        Ok(self.mateid(utils::info_tag_mateid(record)?.map(|e| e.to_vec())))
    }
}

impl Call {
    pub(crate) fn register_heuristic(&mut self, heuristic: Hint) {
        self.hints.insert(heuristic);
    }

    pub(crate) fn write_preprocessed_record(&self, bcf_writer: &mut bcf::Writer) -> Result<()> {
        let rid = bcf_writer.header().name2rid(&self.chrom)?;

        let variant = self.variant.as_ref().unwrap();

        let mut record = bcf_writer.empty_record();
        record.set_rid(Some(rid));
        record.set_pos(self.pos as i64);
        // set ID if present
        if let Some(ref id) = self.id {
            record.set_id(id)?;
        }

        // set alleles
        record.set_alleles(&[&variant.ref_allele[..], &variant.alt_allele[..]])?;

        // set tags
        if let Some(svlen) = variant.svlen {
            record.push_info_integer(b"SVLEN", &[svlen])?;
        }
        if let Some(end) = variant.end {
            record.push_info_integer(b"END", &[end as i32])?;
        }
        if let Some(ref event) = variant.event {
            record.push_info_string(b"EVENT", &[event])?;
        }
        if let Some(ref svtype) = variant.svtype {
            record.push_info_string(b"SVTYPE", &[svtype])?;
        }
        if let Some(ref mateid) = self.mateid {
            record.push_info_string(b"MATEID", &[mateid])?;
        }
        if let Some(heterozygosity) = self.heterozygosity {
            record.push_info_float(
                b"HETEROZYGOSITY",
                &[*PHREDProb::from(heterozygosity) as f32],
            )?;
        }
        if let Some(somatic_effective_mutation_rate) = self.somatic_effective_mutation_rate {
            record.push_info_float(
                b"SOMATIC_EFFECTIVE_MUTATION_RATE",
                &[*PHREDProb::from(somatic_effective_mutation_rate) as f32],
            )?;
        }

        self.write_record_aux_info(variant, &mut record)?;

        self.aux_info.write(&mut record, &OMIT_AUX_INFO)?;

        // set qual
        record.set_qual(f32::missing());

        // add raw observations
        if let Some(ref pileup) = variant.pileup {
            write_observations(pileup, &mut record)?;
        }

        bcf_writer.write(&record)?;

        Ok(())
    }

    fn write_record_aux_info(&self, variant: &Variant, record: &mut bcf::Record) -> Result<()> {
        if let VariantPrecision::Imprecise {
            ref cistart,
            ref ciend,
        } = variant.precision
        {
            record.push_info_flag(b"IMPRECISE")?;
            let get_data = |ci: &RangeInclusive<u64>| [*ci.start() as i32, *ci.end() as i32];
            record.push_info_integer(b"CIPOS", &get_data(cistart))?;
            if let Some(ciend) = ciend {
                record.push_info_integer(b"CIEND", &get_data(ciend))?;
            }
        }
        Ok(())
    }

    pub(crate) fn write_final_record(&mut self, bcf_writer: &mut bcf::Writer) -> Result<()> {
        let rid = bcf_writer.header().name2rid(&self.chrom)?;

        let variant = self.variant.as_ref().unwrap();

        let ref_allele = &variant.ref_allele;
        let mut record = bcf_writer.empty_record();
        record.set_rid(Some(rid));
        record.set_pos(self.pos as i64);
        // set ID if present
        if let Some(ref id) = self.id {
            record.set_id(id)?;
        }

        let mut event_probs = Vec::new();
        let mut allelefreq_estimates = VecMap::new();
        let mut observations = VecMap::new();
        let mut simple_alt_observations = VecMap::new();
        let mut simple_ref_observations = VecMap::new();
        let mut omitted_observations = VecMap::new();
        let mut vaf_densities = VecMap::new();
        let mut obs_counts = VecMap::new();
        let mut strand_bias = VecMap::new();
        let mut read_orientation_bias = VecMap::new();
        let mut read_position_bias = VecMap::new();
        let mut softclip_bias = VecMap::new();
        let mut homopolymer_error = VecMap::new();
        let mut alt_locus_bias = VecMap::new();
        let mut alleles = Vec::new();
        let mut svlens = Vec::new();
        let mut events = Vec::new();
        let mut svtypes = Vec::new();
        let mut ends = Vec::new();
        alleles.push(&ref_allele[..]);

        alleles.push(&variant.alt_allele[..]);

        for (event, prob) in variant
            .event_probs
            .as_ref()
            .expect("bug: event probs must be set")
        {
            event_probs.push((event, *prob));
        }
        event_probs.sort_unstable_by_key(|(_, prob)| OrderedFloat(-prob.0));

        for (i, sample_info) in variant.sample_info.iter().enumerate() {
            if let Some(ref sample_info) = sample_info {
                strand_bias.insert(
                    i,
                    match sample_info.artifacts.strand_bias() {
                        StrandBias::None { .. } => b'.',
                        StrandBias::Forward => b'+',
                        StrandBias::Reverse => b'-',
                    },
                );
                read_orientation_bias.insert(
                    i,
                    match sample_info.artifacts.read_orientation_bias() {
                        ReadOrientationBias::None => b'.',
                        ReadOrientationBias::F1R2 => b'>',
                        ReadOrientationBias::F2R1 => b'<',
                    },
                );
                read_position_bias.insert(
                    i,
                    match sample_info.artifacts.read_position_bias() {
                        ReadPositionBias::None { .. } => b'.',
                        ReadPositionBias::Some => b'^',
                    },
                );
                softclip_bias.insert(
                    i,
                    match sample_info.artifacts.softclip_bias() {
                        SoftclipBias::None => b'.',
                        SoftclipBias::Some => b'$',
                    },
                );
                homopolymer_error.insert(
                    i,
                    match sample_info.artifacts.homopolymer_error() {
                        HomopolymerError::None { .. } => b'.',
                        HomopolymerError::Some { .. } => b'*',
                    },
                );
                alt_locus_bias.insert(
                    i,
                    match sample_info.artifacts.alt_locus_bias() {
                        AltLocusBias::None => b'.',
                        AltLocusBias::Some => b'*',
                    },
                );

                allelefreq_estimates.insert(i, *sample_info.allelefreq_estimate as f32);
                obs_counts.insert(
                    i,
                    expected_depth(sample_info.pileup.read_observations()) as i32,
                );

                observations.insert(
                    i,
                    utils::generalized_cigar(
                        sample_info.pileup.read_observations().iter().map(|obs| {
                            let score = format!("{}", obs.max_bayes_factor());
                            format!(
                                "{}{}{}{}{}{}{}{}{}",
                                if obs.is_max_mapq {
                                    score.to_ascii_uppercase()
                                } else {
                                    score.to_ascii_lowercase()
                                },
                                if let Some(dist) = obs.third_allele_evidence {
                                    format!("{}", dist)
                                } else {
                                    ".".to_owned()
                                },
                                if obs.paired { 'p' } else { 's' },
                                match obs.alt_locus {
                                    AltLocus::Major => '#',
                                    AltLocus::Some => '*',
                                    AltLocus::None => '.',
                                },
                                match obs.strand {
                                    Strand::Both => '*',
                                    Strand::Reverse => '-',
                                    Strand::Forward => '+',
                                    Strand::None => '.',
                                },
                                match obs.read_orientation {
                                    SequenceReadPairOrientation::F1R2 => '>',
                                    SequenceReadPairOrientation::F2R1 => '<',
                                    SequenceReadPairOrientation::None => '*',
                                    _ => '!',
                                },
                                match obs.read_position {
                                    ReadPosition::Major => '^',
                                    ReadPosition::Some => '*',
                                },
                                if obs.softclipped { '$' } else { '.' },
                                if obs.has_homopolymer_error() {
                                    '*'
                                } else {
                                    '.'
                                },
                            )
                        }),
                        false,
                        |(item, _count)| {
                            if item.starts_with('N') {
                                2
                            } else if item.starts_with('E') {
                                1
                            } else {
                                0
                            }
                        },
                    ),
                );

                let fmt_simple_obs = |alt_allele: bool| {
                    utils::generalized_cigar(
                        sample_info
                            .pileup
                            .read_observations()
                            .iter()
                            .filter_map(|obs| {
                                let bf = if alt_allele {
                                    obs.bayes_factor_alt()
                                } else {
                                    obs.bayes_factor_ref()
                                };
                                let score = bayes_factor_to_letter(bf);
                                let keep = (alt_allele
                                    && obs.prob_alt_orig() > obs.prob_ref_orig())
                                    || (!alt_allele && obs.prob_alt_orig() <= obs.prob_ref_orig());
                                if keep {
                                    Some(format!(
                                        "{}",
                                        if obs.is_max_mapq {
                                            score.to_ascii_uppercase()
                                        } else {
                                            score.to_ascii_lowercase()
                                        }
                                    ))
                                } else {
                                    None
                                }
                            }),
                        false,
                        |(item, _count)| {
                            if item.starts_with('R') {
                                2
                            } else if item.ends_with('E') {
                                1
                            } else {
                                0
                            }
                        },
                    )
                };

                simple_alt_observations.insert(i, fmt_simple_obs(true));
                simple_ref_observations.insert(i, fmt_simple_obs(false));

                omitted_observations.insert(i, sample_info.pileup.n_filtered_out_observations());

                vaf_densities.insert(i, sample_info.vaf_dist.clone());
            }
        }

        if let Some(svlen) = variant.svlen {
            svlens.push(svlen);
        }

        if let Some(ref event) = variant.event {
            events.push(event.as_slice());
        }
        if let Some(ref svtype) = variant.svtype {
            svtypes.push(svtype.as_slice());
        }

        if let Some(end) = variant.end {
            ends.push(end as i32);
        }

        // set alleles
        record.set_alleles(&alleles)?;

        if !svlens.is_empty() {
            record.push_info_integer(b"SVLEN", &svlens)?;
        }
        if !svtypes.is_empty() {
            record.push_info_string(b"SVTYPE", &svtypes)?;
        }
        if !events.is_empty() {
            record.push_info_string(b"EVENT", &events)?;
        }
        if !ends.is_empty() {
            record.push_info_integer(b"END", &ends)?;
        }
        self.write_record_aux_info(variant, &mut record)?;

        if let Some(ref mateid) = self.mateid {
            record.push_info_string(b"MATEID", &[mateid])?;
        }

        let is_missing_data = variant
            .sample_info
            .iter()
            .all(|sample| sample.as_ref().is_none_or(|info| info.pileup.is_empty()));

        if is_missing_data {
            self.hints.insert(Hint::MissingData);
        }

        if !self.hints.is_empty() {
            let hints: Vec<&[u8]> = self
                .hints
                .iter()
                .map(|hint| <&Hint as Into<&'static str>>::into(hint).as_bytes())
                .collect();
            record.push_info_string(b"HINTS", &hints)?;
        }

        self.aux_info.write(&mut record, &OMIT_AUX_INFO)?;

        // set qual
        record.set_qual(f32::missing());

        // set event probabilities
        let mut push_prob =
            |event, prob| record.push_info_float(event_tag_name(event).as_bytes(), &[prob]);
        if is_missing_data {
            // missing data
            for (event, _) in event_probs {
                push_prob(event, f32::missing())?;
            }
        } else {
            assert!(
                !event_probs.iter().any(|(_, prob)| prob.is_nan()),
                "bug: event probability is NaN but not all observations are empty for record at {}:{}",
                str::from_utf8(&self.chrom).unwrap(),
                self.pos + 1,
            );
            for (event, prob) in event_probs {
                let prob = PHREDProb::from(prob).abs() as f32;
                push_prob(event, prob)?;
            }
        }

        // set sample info
        if !is_missing_data {
            let dp = obs_counts.values().cloned().collect_vec();
            record.push_format_integer(b"DP", &dp)?;

            let afs = allelefreq_estimates.values().cloned().collect_vec();
            record.push_format_float(b"AF", &afs)?;

            let saobs = simple_alt_observations
                .values()
                .map(|sample_obs| {
                    if sample_obs.is_empty() {
                        b"."
                    } else {
                        sample_obs.as_bytes()
                    }
                })
                .collect_vec();
            record.push_format_string(b"SAOBS", &saobs)?;

            let srobs = simple_ref_observations
                .values()
                .map(|sample_obs| {
                    if sample_obs.is_empty() {
                        b"."
                    } else {
                        sample_obs.as_bytes()
                    }
                })
                .collect_vec();
            record.push_format_string(b"SROBS", &srobs)?;

            let obs = observations
                .values()
                .map(|sample_obs| {
                    if sample_obs.is_empty() {
                        b"."
                    } else {
                        sample_obs.as_bytes()
                    }
                })
                .collect_vec();
            record.push_format_string(b"OBS", &obs)?;

            let oobs = omitted_observations
                .values()
                .map(|n| **n as i32)
                .collect_vec();
            record.push_format_integer(b"OOBS", &oobs)?;

            let sb = strand_bias.values().map(|sb| vec![*sb]).collect_vec();
            record.push_format_string(b"SB", &sb)?;

            let rob = read_orientation_bias
                .values()
                .map(|rob| vec![*rob])
                .collect_vec();
            record.push_format_string(b"ROB", &rob)?;

            let rpb = read_position_bias
                .values()
                .map(|rpb| vec![*rpb])
                .collect_vec();
            record.push_format_string(b"RPB", &rpb)?;

            let scb = softclip_bias.values().map(|scb| vec![*scb]).collect_vec();
            record.push_format_string(b"SCB", &scb)?;

            let he = homopolymer_error.values().map(|he| vec![*he]).collect_vec();
            record.push_format_string(b"HE", &he)?;

            let alb = alt_locus_bias.values().map(|alb| vec![*alb]).collect_vec();
            record.push_format_string(b"ALB", &alb)?;

            let vaf_densities = vaf_densities
                .values()
                .map(|vaf_dist| {
                    vaf_dist.as_ref().map_or_else(
                        || b".".to_vec(),
                        |dist| {
                            dist.iter()
                                .sorted_by_key(|(vaf, _)| *vaf)
                                .map(|(vaf, prob)| {
                                    format!("{:.3}={:.2}", **vaf, *PHREDProb::from(*prob))
                                })
                                .join(",")
                                .into_bytes()
                        },
                    )
                })
                .collect_vec();
            record.push_format_string(b"AFD", &vaf_densities)?;
        } else {
            record.push_format_integer(b"DP", &vec![i32::missing(); variant.sample_info.len()])?;
            record.push_format_float(b"AF", &vec![f32::missing(); variant.sample_info.len()])?;
            record.push_format_string(b"SAOBS", &vec![b".".to_vec(); variant.sample_info.len()])?;
            record.push_format_string(b"SROBS", &vec![b".".to_vec(); variant.sample_info.len()])?;
            record.push_format_string(b"OBS", &vec![b".".to_vec(); variant.sample_info.len()])?;
            record
                .push_format_integer(b"OOBS", &vec![i32::missing(); variant.sample_info.len()])?;
            record.push_format_string(b"SB", &vec![b".".to_vec(); variant.sample_info.len()])?;
            record.push_format_string(b"ROB", &vec![b".".to_vec(); variant.sample_info.len()])?;
            record.push_format_string(b"RPB", &vec![b".".to_vec(); variant.sample_info.len()])?;
            record.push_format_string(b"AFD", &vec![b".".to_vec(); variant.sample_info.len()])?;
        }

        bcf_writer.write(&record)?;
        Ok(())
    }
}

#[derive(Default, Clone, Debug, Builder, Getters)]
pub(crate) struct Variant {
    #[builder(private)]
    ref_allele: Vec<u8>,
    #[builder(private)]
    alt_allele: Vec<u8>,
    #[builder(private, default = "None")]
    svlen: Option<i32>,
    #[builder(private, default = "None")]
    svtype: Option<Vec<u8>>,
    #[builder(private, default = "None")]
    event: Option<Vec<u8>>,
    #[builder(private, default = "None")]
    end: Option<u64>,
    #[builder(private, default)]
    precision: VariantPrecision,
    #[builder(private, default = "None")]
    #[getset(get = "pub(crate)")]
    event_probs: Option<HashMap<String, LogProb>>,
    #[builder(default = "None")]
    pileup: Option<Rc<Pileup>>,
    #[builder(default)]
    #[getset(get = "pub(crate)")]
    sample_info: Vec<Option<SampleInfo>>,
}

impl VariantBuilder {
    /// Build the variant from a single-allele bcf record.
    pub(crate) fn record(&mut self, record: &mut bcf::Record) -> Result<&mut Self> {
        let alleles = record.alleles();
        Ok(self
            .ref_allele(alleles[0].to_owned())
            .alt_allele(alleles[1].to_owned())
            .svlen(record.info(b"SVLEN").integer()?.map(|v| v[0]))
            .event(utils::info_tag_event(record)?.map(|e| e.to_vec()))
            .svtype(utils::info_tag_svtype(record)?.map(|s| s.to_vec()))
            .end(record.info(b"END").integer()?.map(|v| v[0] as u64))
            .precision(VariantPrecision::try_from(&*record)?))
    }

    pub(crate) fn variant(
        &mut self,
        variant: &model::Variant,
        haplotype: &Option<HaplotypeIdentifier>,
        start: usize,
        chrom_seq: Option<&[u8]>,
    ) -> &mut Self {
        self.event(haplotype.as_ref().map(|haplotype| match haplotype {
            HaplotypeIdentifier::Event(event) => event.clone(),
        }));

        match variant {
            model::Variant::Deletion(l) => {
                let l = *l;
                let svlen = -(l as i32);
                if l <= 50 {
                    self.ref_allele(
                        chrom_seq.unwrap()[start..start + 1 + l as usize].to_ascii_uppercase(),
                    )
                    .alt_allele(chrom_seq.unwrap()[start..start + 1].to_ascii_uppercase())
                    .svlen(Some(svlen))
                } else {
                    self.ref_allele(chrom_seq.unwrap()[start..start + 1].to_ascii_uppercase())
                        .alt_allele(b"<DEL>".to_ascii_uppercase())
                        .svlen(Some(svlen))
                        .svtype(Some(b"DEL".to_vec()))
                }
            }
            model::Variant::Insertion(ref seq) => {
                let svlen = seq.len() as i32;
                let ref_allele = vec![chrom_seq.unwrap()[start]];
                let mut alt_allele = ref_allele.clone();
                alt_allele.extend(seq);

                self.ref_allele(ref_allele.to_ascii_uppercase())
                    .alt_allele(alt_allele.to_ascii_uppercase())
                    .svlen(Some(svlen))
                    .svtype(Some(b"INS".to_vec()))
            }
            model::Variant::Snv(base) => self
                .ref_allele(chrom_seq.unwrap()[start..start + 1].to_ascii_uppercase())
                .alt_allele([*base].to_ascii_uppercase()),
            model::Variant::Mnv(bases) => self
                .ref_allele(chrom_seq.unwrap()[start..start + bases.len()].to_ascii_uppercase())
                .alt_allele(bases.to_ascii_uppercase()),
            model::Variant::Breakend {
                ref_allele,
                spec,
                precision,
            } => self
                .ref_allele(ref_allele.to_ascii_uppercase())
                .alt_allele(spec.to_vec())
                .svtype(Some(b"BND".to_vec()))
                .precision(precision.clone()),
            model::Variant::Inversion(len) => self
                .ref_allele(chrom_seq.unwrap()[start..start + 1].to_ascii_uppercase())
                .alt_allele(b"<INV>".to_vec())
                .svtype(Some(b"INV".to_vec()))
                .end(Some(start as u64 + len)), // end tag is inclusive but one-based (hence - 1 + 1)
            model::Variant::Duplication(len) => self
                .ref_allele(chrom_seq.unwrap()[start..start + 1].to_ascii_uppercase())
                .alt_allele(b"<DUP>".to_vec())
                .svtype(Some(b"DUP".to_vec()))
                .end(Some(start as u64 + len)), // end tag is inclusive but one-based (hence - 1 + 1)
            model::Variant::Replacement {
                ref ref_allele,
                ref alt_allele,
            } => self
                .ref_allele(ref_allele.to_owned())
                .alt_allele(alt_allele.to_owned()),
            model::Variant::None => self
                .ref_allele(chrom_seq.unwrap()[start..start + 1].to_ascii_uppercase())
                .alt_allele(b"<REF>".to_ascii_uppercase()),
        }
    }
}

#[derive(Debug, Clone, Builder, Getters, CopyGetters)]
pub(crate) struct SampleInfo {
    #[getset(get_copy = "pub(crate)")]
    allelefreq_estimate: AlleleFreq,
    pileup: Rc<Pileup>,
    artifacts: Artifacts,
    #[getset(get = "pub(crate)")]
    vaf_dist: Option<HashMap<AlleleFreq, LogProb>>,
}

pub(crate) fn chrom<'a>(inbcf: &'a bcf::Reader, record: &bcf::Record) -> &'a [u8] {
    inbcf.header().rid2name(record.rid().unwrap()).unwrap()
}

pub(crate) fn event_tag_name(event: &str) -> String {
    format!("PROB_{}", event.to_ascii_uppercase())
}

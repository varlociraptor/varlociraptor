use anyhow::Result;
use bio_types::genome::{AbstractLocus, Locus};
use progress_logger::ProgressLogger;
use rust_htslib::bcf::{self, Read};
use vec_map::VecMap;

use crate::variants::model;
use crate::{errors, utils};

pub(crate) struct VariantBuffer {
    reader: bcf::Reader,
    variants: Vec<model::Variant>,
    variant_index: usize,
    record_infos: VecMap<RecordInfo>,
    locus: Option<Locus>,
    current_record: Option<bcf::Record>,
    progress_logger: ProgressLogger,
    skips: utils::SimpleCounter<utils::collect_variants::SkipReason>,
    log_each_record: bool,
    record_index: isize,
    state: State,
}

impl VariantBuffer {
    pub(crate) fn new(
        reader: bcf::Reader,
        progress_logger: ProgressLogger,
        log_each_record: bool,
    ) -> Self {
        VariantBuffer {
            reader,
            variants: Vec::new(),
            variant_index: 0,
            record_infos: VecMap::default(),
            locus: None,
            current_record: None,
            progress_logger,
            skips: utils::SimpleCounter::default(),
            log_each_record,
            record_index: -1,
            state: State::Init,
        }
    }

    pub(crate) fn next(&mut self) -> Result<Option<Variants>> {
        if self.skips.total_count() > 0 && self.skips.total_count() % 100 == 0 {
            self.display_skips();
        }
        loop {
            match self.state {
                State::Init => {
                    self.next_record()?;
                    if self.current_record.is_none() {
                        return Ok(None); // empty file
                    } else {
                        self.state = State::InitLocus;
                    }
                }
                State::InitLocus => {
                    self.variants.clear();
                    self.variant_index = 0;
                    self.locus = Some(self.locus(self.current_record.as_ref().unwrap()));
                    self.add_variants()?;
                    self.state = State::LocusInProgress;
                }
                State::LocusInProgress => {
                    self.next_record()?;
                    if self.current_record.is_none() {
                        // EOF
                        self.state = State::LocusComplete;
                    } else {
                        let previous_locus = self.locus.as_ref().unwrap();
                        let locus = self.locus(self.current_record.as_ref().unwrap());
                        if locus != *previous_locus {
                            if locus.contig() == previous_locus.contig()
                                && locus.pos() < previous_locus.pos()
                            {
                                // unsorted input file, fail with an error
                                return Err(errors::Error::UnsortedVariantFile.into());
                            }
                            self.state = State::LocusComplete;
                        } else {
                            self.add_variants()?;
                        }
                    }
                }
                State::LocusComplete => {
                    if self.variant_index < self.variants.len() {
                        // yield next variant
                        let variant_index = self.variant_index;
                        self.variant_index += 1;

                        let variants = Variants {
                            variant_of_interest: &self.variants[variant_index],
                            before: &self.variants[..variant_index],
                            after: if variant_index < self.variants.len() {
                                &self.variants[variant_index + 1..]
                            } else {
                                &self.variants[..0] // generate empty slice
                            },
                            locus: self.locus.clone().unwrap(),
                            record_info: self.record_infos.get(variant_index).unwrap(),
                        };
                        if self.log_each_record {
                            info!(
                                "Alt variants at previous locus: {}",
                                variants.n_alt_variants()
                            );
                        }

                        return Ok(Some(variants));
                    } else if self.current_record.is_some() {
                        self.state = State::InitLocus;
                    } else {
                        self.display_skips();
                        return Ok(None); // end of file
                    }
                }
            }
        }
    }

    fn next_record(&mut self) -> Result<()> {
        let mut record = self.reader.empty_record();
        match self.reader.read(&mut record) {
            None => {
                self.current_record = None;
            }
            Some(res) => {
                res?;
                self.record_index += 1;
                if self.log_each_record {
                    info!(
                        "Processing record {} at {}:{}",
                        self.record_index,
                        record.contig(),
                        record.pos() + 1,
                    );
                }
                self.progress_logger.update(1u64);
                self.current_record = Some(record);
            }
        }
        Ok(())
    }

    fn add_variants(&mut self) -> Result<()> {
        let record = self.current_record.as_mut().unwrap();
        let variants = utils::collect_variants(record, true, Some(&mut self.skips))?;
        let record_info = RecordInfo::new(
            self.record_index as usize,
            record.id(),
            utils::info_tag_mateid(record).unwrap_or(None),
        );
        for index in self.variants.len()..self.variants.len() + variants.len() {
            self.record_infos.insert(index, record_info.clone());
        }

        self.variants.extend(variants);

        Ok(())
    }

    fn display_skips(&self) {
        for (reason, &count) in self.skips.iter() {
            if count > 0 && count % 100 == 0 {
                info!("Skipped {} {}.", count, reason);
            }
        }
    }

    fn locus(&self, record: &bcf::Record) -> Locus {
        Locus::new(record.contig().to_owned(), record.pos() as u64)
    }
}

#[derive(CopyGetters, Getters)]
pub(crate) struct Variants<'a> {
    #[getset(get_copy = "pub(crate)")]
    variant_of_interest: &'a model::Variant,
    before: &'a [model::Variant],
    after: &'a [model::Variant],
    #[getset(get = "pub(crate)")]
    locus: Locus,
    #[getset(get_copy = "pub(crate)")]
    record_info: &'a RecordInfo,
}

impl<'a> Variants<'a> {
    pub(crate) fn alt_variants(&self) -> impl Iterator<Item = &'a model::Variant> {
        self.before.iter().chain(self.after.iter())
    }

    pub(crate) fn n_alt_variants(&self) -> usize {
        self.before.len() + self.after.len()
    }
}

#[derive(Clone, Debug, new, Getters, CopyGetters)]
pub(crate) struct RecordInfo {
    #[getset(get_copy = "pub(crate)")]
    index: usize,
    #[getset(get = "pub(crate)")]
    id: Vec<u8>,
    #[getset(get = "pub(crate)")]
    mateid: Option<Vec<u8>>,
}

#[derive(Debug, Clone, Copy)]
enum State {
    Init,
    InitLocus,
    LocusInProgress,
    LocusComplete,
}

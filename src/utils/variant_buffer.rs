use std::cell::Cell;

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
    next_locus_record: Option<Cell<bcf::Record>>,
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
            next_locus_record: None,
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
            dbg!(self.state);
            match self.state {
                State::Init => {
                    self.next_locus_record = self.next_record()?;
                    if self.next_locus_record.is_none() {
                        return Ok(None); // empty file
                    } else {
                        self.state = State::InitLocus;
                    }
                }
                State::InitLocus => {
                    self.variants.clear();
                    self.variant_index = 0;
                    self.locus = Some(self.locus(self.next_locus_record.as_ref().unwrap()));
                    self.add_variants(self.next_locus_record.as_mut().unwrap())?;
                    self.state = State::LocusInProgress;
                }
                State::LocusInProgress => {
                    if let Some(mut record) = self.next_record()? {
                        let locus = self.locus(&record);

                        if locus != *self.locus.as_ref().unwrap() {
                            self.state = State::LocusComplete;
                            self.next_locus_record = Some(record);
                        } else {
                            self.add_variants(&mut record)?;
                        }
                    } else {
                        self.state = State::LocusComplete;
                    }
                }
                State::LocusComplete => {
                    if self.variant_index < self.variants.len() {
                        // yield next variant
                        let variant_index = self.variant_index;
                        self.variant_index += 1;

                        return Ok(Some(Variants {
                            variant_of_interest: &self.variants[variant_index],
                            before: &self.variants[..variant_index],
                            after: if variant_index < self.variants.len() {
                                &self.variants[variant_index + 1..]
                            } else {
                                &self.variants[..0] // generate empty slice
                            },
                            locus: self.locus.clone().unwrap(),
                            record_info: self.record_infos.get(variant_index).unwrap(),
                        }));
                    } else if self.next_locus_record.is_some() {
                        self.state = State::InitLocus;
                    } else {
                        self.display_skips();
                        return Ok(None); // end of file
                    }
                }
            }
        }
    }
    
    fn next_record(&mut self) -> Result<Option<bcf::Record>> {
        let mut record = self.reader.empty_record();
        match self.reader.read(&mut record) {
            None => {
                Ok(None)
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
                Ok(Some(record))
            },
        }
    }

    // pub(crate) fn next(&mut self) -> Result<Option<Variants>> {
    //     loop {
    //         if self.variants.is_empty() {
    //             // Step 1: start new locus
    //             if let Some(mut record) = self.next_locus_record.take() {
    //                 self.add_variants(&mut record)?;
    //                 self.locus = Some(self.locus(&record));
    //                 self.next_locus_record = None;
    //             }
    //         }

    //         let mut end_of_file = false;

    //         // Step 2: read until next locus is reached
    //         loop {
    //             let mut record = self.reader.empty_record();
    //             match self.reader.read(&mut record) {
    //                 None => {
    //                     end_of_file = true;
    //                     break;
    //                 }
    //                 Some(res) => res?,
    //             }
    //             let record_locus = self.locus(&record);
    //             self.progress_logger.update(1u64);
    //             dbg!(&self.variants);

    //             if self.log_each_record {
    //                 info!(
    //                     "Processing record at {}:{}",
    //                     record.contig(),
    //                     record.pos() + 1
    //                 );
    //             }

    //             if let Some(ref current_locus) = self.locus {
    //                 if record_locus != *current_locus {
    //                     if record_locus.contig() == current_locus.contig()
    //                         && record_locus.pos() > current_locus.pos()
    //                     {
    //                         return Err(errors::Error::UnsortedVariantFile.into());
    //                     }
    //                     // new locus, store and stop
    //                     self.next_locus_record = Some(record);
    //                     self.record_index += 1;
    //                     break;
    //                 }
    //             } else {
    //                 // no locus yet, record this one
    //                 self.locus = Some(record_locus);
    //             }

    //             self.add_variants(&mut record)?;

    //             if self.skips.total_count() > 0 && self.skips.total_count() % 100 == 0 {
    //                 self.display_skips();
    //             }
    //             self.record_index += 1;
    //         }
    //         dbg!(&self.variants);

    //         if let Some(ref locus) = self.locus {
    //             // Step 3: return slices
    //             let variant_index = self.variant_index;
    //             dbg!((variant_index, self.variants.len()));
    //             self.variant_index += 1; // move to next for next call

    //             if variant_index >= self.variants.len() {
    //                 // all variants returned, clear and prepare for next locus
    //                 self.variants.clear();

    //                 if variant_index > 0 {
    //                     self.variant_index = 0;
                        
    //                 }
    //                 continue;
    //             }

    //             if end_of_file {
    //                 self.display_skips();
    //             }
    //             dbg!(&self.variants);

    //             return Ok(Some(Variants {
    //                 variant_of_interest: &self.variants[variant_index],
    //                 before: &self.variants[..variant_index],
    //                 after: if variant_index < self.variants.len() {
    //                     &self.variants[variant_index + 1..]
    //                 } else {
    //                     &self.variants[..0] // generate empty slice
    //                 },
    //                 locus: locus.clone(),
    //                 record_info: self.record_infos.get(variant_index).unwrap(),
    //             }));
    //         } else {
    //             if end_of_file {
    //                 self.display_skips();
    //             }
    //             return Ok(None);
    //         }
    //     }
    // }

    fn add_variants(&mut self, record: &mut bcf::Record) -> Result<()> {
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
    #[getset(get_copy = "pub(crate)")]
    before: &'a [model::Variant],
    #[getset(get_copy = "pub(crate)")]
    after: &'a [model::Variant],
    #[getset(get = "pub(crate)")]
    locus: Locus,
    #[getset(get_copy = "pub(crate)")]
    record_info: &'a RecordInfo,
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
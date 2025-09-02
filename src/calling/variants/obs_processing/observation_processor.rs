use rust_htslib::bcf;

use crate::{
    calling::variants::Call,
    utils::{aux_info::AuxInfoCollector, variant_buffer::Variants},
    variants::evidence::observations::pileup::Pileup,
    variants::sample::Sample,
};
use anyhow::Result;

pub trait ObservationProcessor {
    type Realigner;

    fn writer(&self, aux_info_collector: &AuxInfoCollector) -> Result<bcf::Writer>;

    fn process(&mut self) -> Result<()>;

    fn write_observations(&self, pileup: &Pileup, variants: &Variants) -> Result<()>;

    fn process_variant(&self, variants: Variants, sample: &mut Sample) -> Result<Vec<Call>>;

    fn process_pileup(&self, variants: &Variants, sample: &mut Sample) -> Result<Option<Pileup>>;
}

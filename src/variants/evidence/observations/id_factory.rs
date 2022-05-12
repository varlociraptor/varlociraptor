use std::collections::HashMap;

use bio_types::genome::AbstractLocus;

use super::read_observation::Evidence;

#[derive(Default, Debug)]
pub(crate) struct ObservationIdFactory {
    last_ends: HashMap<String, u64>,
    ids: HashMap<Vec<u8>, u64>,
    next_id: u64,
}

impl ObservationIdFactory {
    pub(crate) fn register<E>(&mut self, evidence: &E) -> u64
    where
        E: Evidence,
    {
        // check whether ids can be reset
        let start_locus = evidence.start_locus();
        if let Some(last_end) = self.last_ends.get(start_locus.contig()) {
            if *last_end < start_locus.pos() {
                self.last_ends.remove(start_locus.contig());
            }
        }
        if self.last_ends.is_empty() {
            // no overlapping end remaining, set next_id to 0
            self.next_id = 0;
        }

        // record end
        let end_locus = evidence.end_locus();
        let last_end = self
            .last_ends
            .entry(end_locus.contig().to_owned())
            .or_default();
        if *last_end < end_locus.pos() {
            *last_end = end_locus.pos();
        }

        let name = evidence.name();

        if let Some(id) = self.ids.get(name) {
            *id
        } else {
            let id = self.next_id;
            self.ids.insert(name.to_owned(), id);
            id
        }
    }
}

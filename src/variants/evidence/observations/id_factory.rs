use std::collections::HashMap;

use bio_types::genome::AbstractLocus;

use super::read_observation::Evidence;

#[derive(Default, Debug)]
pub(crate) struct ObservationIdFactory {
    ids: HashMap<Vec<u8>, u64>,
    next_id: u64,
    current_contig: String,
}

impl ObservationIdFactory {
    pub(crate) fn register_contig(&mut self, contig: &str) {
        if self.current_contig != contig {
            self.ids.clear();
            self.next_id = 0;
            self.current_contig = contig.to_owned();
        }
    }
    pub(crate) fn register<E>(&mut self, evidence: &E) -> u64
    where
        E: Evidence,
    {
        if self.ids.contains_key(evidence.name()) {
            *self.ids.get(evidence.name()).unwrap()
        } else {
            let id = self.next_id;
            self.ids.insert(evidence.name().to_vec(), id);
            self.next_id += 1;
            id
        }
    }
}

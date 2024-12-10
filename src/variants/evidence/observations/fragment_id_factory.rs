use std::collections::HashMap;

use super::read_observation::{Evidence, EvidenceIdentifier};

#[derive(Default, Debug)]
pub(crate) struct FragmentIdFactory {
    ids: HashMap<EvidenceIdentifier, u64>,
    next_id: u64,
    current_contig: String,
}

impl FragmentIdFactory {
    pub(crate) fn register_contig(&mut self, contig: &str) {
        if self.current_contig != contig {
            self.ids.clear();
            self.next_id = 0;
            self.current_contig = contig.to_owned();
        }
    }
    pub(crate) fn register(&mut self, evidence: &Evidence) -> u64 {
        if self.ids.contains_key(&evidence.id()) {
            *self.ids.get(&evidence.id()).unwrap()
        } else {
            let id = self.next_id;
            self.ids.insert(evidence.id(), id);
            self.next_id += 1;
            id
        }
    }
}

use std::num::ParseIntError;

use bio::io::fasta::Sequence;

#[derive(Clone, Debug)]
pub struct Batch {
    idx: usize,
    total: usize
}

impl From<&str> for Batch {
    pub fn from(v: &str) -> Result<Self, ParseIntError> {
        let (idx, total) = v.split_at("/");
        Batch {
            idx: idx.parse()?,
            total: total.parse()?,
        }
    }
}

impl Batch {
    pub fn is_valid(&self) -> bool {
        self.idx < self.total
    }
}

/// Select sequences for a given batch. We do this by solving the multi-way partitioning problem
/// with the well-known naive greedy heuristic.
pub fn sequence_batch<'a>(batch: &Batch, mut sequences: Vec<Sequence>) -> Vec<Sequence> {
    sequences.sort_by_key(|seq| - seq.len as i64);
    let mut k = 0;
    let mut batch_seqs = Vec::new();
    for seq in sequences {
        if k == batch.idx {
            // add sequence to batch of interest
            batch_seqs.push(seq);
        }
        k += 1;
        if k >= batch.total {
            // start with first batch again
            k = 0;
        }
    }

    batch_seqs
}

use std::fs;
use std::str;
use std::sync::Arc;
use std::sync::{Mutex, RwLock};

use anyhow::Result;
use bio::io::fasta;
use lru_time_cache::LruCache;

/// A lazy buffer for reference sequences.
pub(crate) struct Buffer {
    reader: RwLock<fasta::IndexedReader<fs::File>>,
    sequences: Mutex<LruCache<String, Arc<Vec<u8>>>>,
}

impl Buffer {
    pub(crate) fn new(fasta: fasta::IndexedReader<fs::File>, capacity: usize) -> Self {
        Buffer {
            reader: RwLock::new(fasta),
            sequences: Mutex::new(LruCache::with_capacity(capacity)),
        }
    }

    pub(crate) fn sequences(&self) -> Vec<fasta::Sequence> {
        self.reader.read().unwrap().index.sequences()
    }

    /// Load given chromosome and return it as a slice. This is O(1) if chromosome was loaded before.
    pub(crate) fn seq(&self, chrom: &str) -> Result<Arc<Vec<u8>>> {
        let mut sequences = self.sequences.lock().unwrap();

        if !sequences.contains_key(chrom) {
            let mut sequence = Arc::new(Vec::new());
            {
                let mut reader = self.reader.write().unwrap();
                reader.fetch_all(chrom)?;
                reader.read(Arc::get_mut(&mut sequence).unwrap())?;
            }

            sequences.insert(chrom.to_owned(), Arc::clone(&sequence));
            Ok(sequence)
        } else {
            Ok(Arc::clone(sequences.get(chrom).unwrap()))
        }
    }
}

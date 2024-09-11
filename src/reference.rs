use std::fs;
use std::path::{Path, PathBuf};
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
    reference_path: Option<PathBuf>,
}

impl Buffer {
    pub(crate) fn from_path<P: AsRef<Path> + std::fmt::Debug>(
        path: P,
        capacity: usize,
    ) -> Result<Self> {
        let fasta: fasta::IndexedReader<fs::File> = fasta::IndexedReader::from_file(&path)?;
        Ok(Buffer {
            reader: RwLock::new(fasta),
            sequences: Mutex::new(LruCache::with_capacity(capacity)),
            reference_path: Some(path.as_ref().to_path_buf()),
        })
    }

    pub(crate) fn reference_path(&self) -> Option<&PathBuf> {
        self.reference_path.as_ref()
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

use std::fs;
use std::str;
use std::sync::Arc;

use anyhow::Result;
use bio::io::fasta;

/// A lazy buffer for reference sequences.
pub struct Buffer {
    pub(crate) reader: fasta::IndexedReader<fs::File>,
    chrom: Option<Vec<u8>>,
    sequence: Arc<Vec<u8>>,
}

impl Buffer {
    pub fn new(fasta: fasta::IndexedReader<fs::File>) -> Self {
        Buffer {
            reader: fasta,
            chrom: None,
            sequence: Arc::new(Vec::new()),
        }
    }

    fn is_current_chrom(&self, chrom: &[u8]) -> bool {
        if let Some(ref last_chrom) = self.chrom {
            if last_chrom == &chrom {
                return true;
            }
        }
        false
    }

    /// Load given chromosome and return it as a slice. This is O(1) if chromosome was loaded before.
    pub fn seq(&mut self, chrom: &[u8]) -> Result<Arc<Vec<u8>>> {
        if self.is_current_chrom(chrom) {
            return Ok(Arc::clone(&self.sequence));
        }

        self.reader.fetch_all(str::from_utf8(chrom)?)?;
        self.reader
            .read(Arc::get_mut(&mut self.sequence).expect(
                "bug: reference on sequence buffer still alive while trying to refill it",
            ))?;

        self.chrom = Some(chrom.to_owned());

        Ok(Arc::clone(&self.sequence))
    }
}

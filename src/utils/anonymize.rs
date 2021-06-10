use std::collections::HashMap;

use anyhow::Result;
use rand::{seq::SliceRandom, thread_rng};
use rust_htslib::{bam, bcf};
use uuid::Uuid;

pub(crate) struct Anonymizer {
    replacements: HashMap<u8, u8>,
    qnames: HashMap<Vec<u8>, Vec<u8>>,
}

impl Anonymizer {
    pub(crate) fn new() -> Self {
        let mut alphabet = b"ACGT".to_owned();
        alphabet.shuffle(&mut thread_rng());
        let mut replacements = HashMap::new();
        for (base, repl) in b"ACGT".into_iter().zip(&alphabet) {
            replacements.insert(*base, *repl);
        }
        replacements.insert(b'N', b'N');

        Anonymizer {
            replacements,
            qnames: HashMap::new(),
        }
    }

    pub(crate) fn anonymize_seq(&self, seq: &mut [u8]) {
        for base in seq.iter_mut() {
            *base = self.replacements[base];
        }
    }

    pub(crate) fn anonymize_bcf_record(&self, record: &mut bcf::Record) -> Result<()> {
        let rand_alleles: Vec<Vec<u8>> = record
            .alleles()
            .iter()
            .map(|allele| {
                if allele
                    .iter()
                    .all(|base| self.replacements.contains_key(base))
                {
                    // replace alleles with chiffre
                    allele.iter().map(|base| self.replacements[base]).collect()
                } else {
                    // non ACGTN allele, keep
                    allele.to_vec()
                }
            })
            .collect();
        let refs: Vec<_> = rand_alleles
            .iter()
            .map(|allele| allele.as_slice())
            .collect();
        record.set_alleles(refs.as_slice())?;

        Ok(())
    }

    pub(crate) fn anonymize_bam_record(&mut self, record: &mut bam::Record) {
        // replace seq with chiffre
        let mut seq = record.seq().as_bytes();
        self.anonymize_seq(&mut seq);
        let qual = record.qual().to_vec();

        // replace qname with uuid, keep pairs in sync by memoizing qnames
        let qname = self
            .qnames
            .entry(record.qname().to_vec())
            .or_insert_with(|| {
                Uuid::new_v4()
                    .to_hyphenated()
                    .to_string()
                    .as_bytes()
                    .to_owned()
            });

        // store in record
        record.set(qname, Some(&record.cigar().take()), &seq, &qual);
    }
}

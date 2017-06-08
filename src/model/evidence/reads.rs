use pairhmm;


/// Calculate read evindence for an indel.
pub struct IndelEvidence {
    gap_params: IndelGapParams,
    pairhmm: PairHMM
}


impl IndelEvidence {
    pub fn new(
        prob_insertion_artifact: LogProb,
        prob_deletion_artifact: LogProb,
        prob_insertion_extend_artifact: LogProb,
        prob_deletion_extend_artifact: LogProb
    ) -> Self {
        IndelEvidence {
            gap_params: IndelGapParams {
                prob_insertion_artifact: prob_insertion_artifact,
                prob_deletion_artifact: prob_deletion_artifact,
                prob_insertion_extend_artifact: prob_insertion_extend_artifact,
                prob_deletion_extend_artifact: prob_deletion_extend_artifact
            },
            pairhmm: pairhmm:PairHMM::new()
        }
    }

    pub fn prob(&mut self,
        record: &bam::Record,
        cigar: &CigarStringView,
        start: u32,
        variant: &Variant,
        ref_seq: &[u8]
    ) -> Result<(LogProb, LogProb), Box<Error>> {
        let read_seq = record.seq();
        let read_qual = record.qual();

        let (read_offset, read_end, breakpoint) = {
            let (varstart, varend) = match variant {
                &Variant::Deletion(_) => (start, start + variant.len()),
                &Variant::Insertion(_) => (start, start + 1),
                &Variant::SNV(_) => panic!("bug: unsupported variant")
            };

            match (
                cigar.read_pos(varstart, true, true)?,
                cigar.read_pos(varend, true, true)?
            ) {
                // read encloses variant
                (Some(qstart), Some(qend)) => {
                    let qstart = qstart as usize;
                    let qend = qend as usize;
                    let read_offset = qstart.saturating_sub(self.indel_haplotype_window as usize);
                    let read_end = cmp::min(
                        qend + self.indel_haplotype_window as usize,
                        read_seq.len()
                    );
                    (read_offset, read_end, varstart as usize)
                },
                (Some(qstart), None) => {
                    let qstart = qstart as usize;
                    let read_offset = qstart.saturating_sub(self.indel_haplotype_window as usize);
                    let read_end = cmp::min(
                        qstart + self.indel_haplotype_window as usize,
                        read_seq.len()
                    );
                    (read_offset, read_end, varstart as usize)
                },
                (None, Some(qend)) => {
                    let qend = qend as usize;
                    let read_offset = qend.saturating_sub(self.indel_haplotype_window as usize);
                    let read_end = cmp::min(
                        qend + self.indel_haplotype_window as usize,
                        read_seq.len()
                    );
                    (read_offset, read_end, varend as usize)
                },
                (None, None) => {
                    panic!(
                        "bug: read does not overlap breakpoint: pos={}, cigar={}, start={}, len={}",
                        record.pos(),
                        cigar,
                        start,
                        variant.len()
                    );
                }
            }
        };

        let start = start as usize;
        // the window on the reference should be a bit larger to allow some flexibility with close
        // indels. But it should not be so large that the read can align outside of the breakpoint.
        let ref_window = (self.indel_haplotype_window as f64 * 1.5) as usize;

        // ref allele
        let prob_ref = self.pairhmm.prob_related(
            &self.gap_params,
            &ReferenceEmissionParams {
                ref_seq: ref_seq,
                read_seq: read_seq,
                read_qual: read_qual,
                read_offset: read_offset,
                read_end: read_end
                ref_offset: breakpoint.saturating_sub(ref_window),
                ref_end: cmp::min(breakpoint + ref_window, ref_seq.len()),
            }
        );

        // alt allele
        let prob_alt = match variant {
            &Variant::Deletion(_) => {
                self.pairhmm.prob_related(
                    &self.gap_params,
                    &DeletionEmissionParams {
                        ref_seq: ref_seq,
                        read_seq: read_seq,
                        read_qual: read_qual,
                        read_offset: read_offset,
                        read_end: read_end,
                        ref_offset: start.saturating_sub(ref_window),
                        ref_end: cmp::min(start + ref_window, ref_seq.len()),
                        del_start: start,
                        del_len: variant.len() as usize
                    }
                )
            },
            &Variant::Insertion(ref ins_seq) => {
                let l = ins_seq.len() as usize;
                self.pairhmm.prob_related(
                    &self.gap_params,
                    &InsertionEmissionParams {
                        ref_seq: ref_seq,
                        read_seq: read_seq,
                        read_qual: read_qual,
                        read_offset: read_offset,
                        read_end: read_end,
                        ref_offset: start.saturating_sub(ref_window),
                        ref_end: cmp::min(start + l + ref_window, ref_seq.len()),
                        ins_start: start,
                        ins_len: l,
                        ins_seq: ins_seq
                    }
                )
            },
            _ => {
                panic!("bug: unsupported variant");
            };
            Ok((prob_ref, prob_alt))
        }
    }
}


/// Calculate probability of read_base given ref_base.
pub fn prob_read_base(read_base: u8, ref_base: u8, base_qual: u8) -> LogProb {
    let prob_miscall = prob_read_base_miscall(base_qual);

    if read_base.to_ascii_uppercase() == ref_base.to_ascii_uppercase() {
        prob_miscall.ln_one_minus_exp()
    } else {
        // TODO replace the second term with technology specific confusion matrix
        prob_miscall + *PROB_CONFUSION
    }
}


/// Calculate probability of read_base given ref_base.
pub fn prob_read_base_miscall(base_qual: u8) -> LogProb {
    LogProb::from(PHREDProb::from((base_qual) as f64))
}


/// Gap parameters for PairHMM.
pub struct IndelGapParams {
    pub prob_insertion_artifact: LogProb,
    pub prob_deletion_artifact: LogProb,
    pub prob_insertion_extend_artifact: LogProb,
    pub prob_deletion_extend_artifact: LogProb
}


impl pairhmm::SemiglobalAlignmentParameters for IndelGapParams {
    #[inline]
    fn prob_gap_x(&self) -> LogProb {
        self.prob_insertion_artifact
    }

    #[inline]
    fn prob_gap_y(&self) -> LogProb {
        self.prob_deletion_artifact
    }

    #[inline]
    fn prob_gap_x_extend(&self) -> LogProb {
        self.prob_insertion_extend_artifact
    }

    #[inline]
    fn prob_gap_y_extend(&self) -> LogProb {
        self.prob_deletion_extend_artifact
    }
}


pub trait EmissionParams: pairhmm::EmissionParameters {
    #[inline]
    fn prob_emit_xy(&self, i: usize, j: usize) -> LogProb {
        let r = self.ref_base(i);
        let j_ = self.project_j(j);
        prob_read_base(self.read_seq[j_], r, self.read_qual[j_])
    }

    #[inline]
    fn prob_emit_x(&self, _: usize) -> LogProb {
        LogProb::ln_one()
    }

    #[inline]
    fn prob_emit_y(&self, j: usize) -> LogProb {
        prob_read_base_miscall(self.read_qual[self.project_j(j)])
    }

    fn ref_base(&self, i: usize) -> u8;

    fn project_j(&self, j: usize) -> usize;
}


/// Emission parameters for PairHMM over reference allele.
pub struct ReferenceEmissionParams<'a> {
    ref_seq: &'a [u8],
    read_seq: &'a [u8],
    read_qual: &'a [u8],
    read_offset: usize,
    ref_offset: usize,
    read_end: usize,
    ref_end: usize
}


impl<'a> EmissionParams for ReferenceEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        self.ref_seq[i + self.ref_offset]
    }

    #[inline]
    fn project_j(&self, j: usize) -> usize {
        j + self.read_offset
    }

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }

    #[inline]
    fn len_y(&self) -> usize {
        self.read_end - self.read_offset
    }
}


/// Emission parameters for PairHMM over deletion allele.
pub struct DeletionEmissionParams<'a> {
    ref_seq: &'a [u8],
    read_seq: &'a [u8],
    read_qual: &'a [u8],
    read_offset: usize,
    ref_offset: usize,
    read_end: usize,
    ref_end: usize,
    del_start: usize,
    del_len: usize
}


impl<'a> EmissionParams for DeletionEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;
        if i_ <= self.del_start {
            self.ref_seq[i_]
        } else {
            self.ref_seq[i_ + self.del_len]
        }
    }

    #[inline]
    fn project_j(&self, j: usize) -> usize {
        j + self.read_offset
    }


    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }

    #[inline]
    fn len_y(&self) -> usize {
        self.read_end - self.read_offset
    }
}


/// Emission parameters for PairHMM over insertion allele.
pub struct InsertionEmissionParams<'a> {
    ref_seq: &'a [u8],
    read_seq: &'a [u8],
    read_qual: &'a [u8],
    read_offset: usize,
    ref_offset: usize,
    read_end: usize,
    ref_end: usize,
    ins_start: usize,
    ins_end: usize,
    ins_len: usize
}


impl<'a> EmissionParams for InsertionEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;
        if i_ <= self.ins_start {
            self.ref_seq[i_]
        } else if i_ > self.ins_end {
            self.ref_seq[i_ - self.ins_len]
        } else {
            self.ins_seq[i_ - (self.ins_start + 1)]
        }
    }

    #[inline]
    fn project_j(&self, j: usize) -> usize {
        j + self.read_offset
    }

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }

    #[inline]
    fn len_y(&self) -> usize {
        self.read_end - self.read_offset
    }
}

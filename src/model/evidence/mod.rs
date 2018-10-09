use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar, CigarStringView};

pub mod fragments;
pub mod observation;
pub mod reads;

pub use self::observation::Evidence;
pub use self::observation::Observation;

pub struct Clips {
    hard: u32,
    soft: u32,
}

impl Clips {
    fn new(cigar: &CigarStringView, trailing: bool) -> Self {
        let ops = if trailing {
            if cigar.len() > 1 {
                (cigar.get(cigar.len() - 1), cigar.get(cigar.len() - 2))
            } else {
                (cigar.get(cigar.len() - 1), None)
            }
        } else {
            (cigar.get(0), cigar.get(1))
        };

        // first comes outer, second comes inner operation
        match ops {
            (Some(&Cigar::HardClip(j)), Some(&Cigar::SoftClip(i))) => Clips { hard: j, soft: i },
            (Some(&Cigar::SoftClip(i)), _) => Clips { hard: 0, soft: i },
            (Some(&Cigar::HardClip(j)), _) => Clips { hard: j, soft: 0 },
            _ => Clips { hard: 0, soft: 0 },
        }
    }

    /// Trailing clips in the cigar string.
    ///
    /// # Arguments
    ///
    /// * `cigar` - the cigar string
    pub fn trailing(cigar: &CigarStringView) -> Self {
        Clips::new(cigar, true)
    }

    /// Leading clips in the cigar string.
    ///
    /// # Arguments
    ///
    /// * `cigar` - the cigar string
    pub fn leading(cigar: &CigarStringView) -> Self {
        Clips::new(cigar, false)
    }

    /// Number of hard clipped bases.
    pub fn hard(&self) -> u32 {
        self.hard
    }

    /// Number of soft clipped bases.
    pub fn soft(&self) -> u32 {
        self.soft
    }

    /// Number of both hard and soft clipped bases.
    pub fn both(&self) -> u32 {
        self.hard + self.soft
    }
}

/// Return maximum indel operation in given cigar string.
pub fn max_indel(cigar: &CigarStringView) -> u32 {
    cigar
        .iter()
        .map(|op| match op {
            &Cigar::Ins(l) => l,
            &Cigar::Del(l) => l,
            _ => 0,
        }).max()
        .unwrap_or(0)
}

/// Calculate the full read length including hard clipped bases.
pub fn read_len(record: &bam::Record, cigar: &CigarStringView) -> u32 {
    record.seq().len() as u32 + Clips::trailing(cigar).hard() + Clips::leading(cigar).hard()
}

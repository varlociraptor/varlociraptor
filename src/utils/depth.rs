use std::collections::BTreeMap;

use bio::stats::{LogProb, Prob};
use rgsl::randist::poisson::{poisson_P, poisson_pdf};
use rust_htslib::bam;

use crate::utils::SimpleCounter;


#[derive(Clone, Debug, Serialize, Deserialize, Default, CopyGetters)]
#[getset(get_copy = "pub(crate)")]
pub(crate) struct Depth {
    mean: f64,
}

impl Depth {
    pub(crate) fn prob_true_pileup(&self, observed_depth: u32) -> LogProb {
        let max_putative_untrue_depth = self.mean.ceil() as u32 - 1;
        if observed_depth > max_putative_untrue_depth {
            LogProb::ln_one()
        } else {
            let total = poisson_P(max_putative_untrue_depth, self.mean);
            LogProb::from(Prob(poisson_pdf(observed_depth, self.mean) / total))
        }
    }
}

#[derive(Clone, Debug, Default)]
pub(crate) struct DepthEstimator {
    buffer: BTreeMap<u64, u64>,
    depths: SimpleCounter<u64>,
    current_tid: Option<i32>,
}

impl DepthEstimator {
    fn update_depths(&mut self) {
        for (pos, depth) in &self.buffer {
            self.depths.incr(*depth);
        }
    }

    /// Takes a mapped record to update the depth estimate.
    pub(crate) fn update(&mut self, record: &bam::Record) {
        assert!(!record.is_unmapped());
        if record.mapq() < 60 {
            return;
        }
        
        // move on to new interval and process anything before
        if let Some(tid) = self.current_tid {
            if tid != record.tid() {
                // new chromosome
                self.update_depths();
                self.buffer.clear();
            } else {
                // new position, same chromosome
                let to_keep = self.buffer.split_off(&(record.pos() as u64));
                self.update_depths();
                self.buffer = to_keep;
            }
        }
        

        // increase depth of covered area
        for pos in record.pos() as u64..record.cigar_cached().expect("bug: cache_cigar() has to be called on the given record").end_pos() as u64 {
            *self.buffer.entry(pos).or_insert(0) += 1;
        }
    }

    /// Estimate depth given the collected records.
    pub(crate) fn estimate(&mut self) -> Depth {
        self.update_depths();
        self.buffer.clear();
        let n = self.depths.values().sum::<usize>() as f64;
        let mean = self.depths.iter().map(|(depth, count)| *depth as usize * *count).sum::<usize>() as f64 / n;

        Depth {
            mean,
        }
    }
}
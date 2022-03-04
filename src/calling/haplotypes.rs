use crate::haplotypes::model::{Data, Likelihood, Marginal, Posterior, Prior};
use crate::variants::model::AlleleFreq;
use anyhow::Result;
use bio::stats::{bayesian::model::Model, probs::LogProb, PHREDProb};
use bv::BitVec;
use derive_builder::Builder;
use hdf5;
use ordered_float::NotNan;
use rust_htslib::bcf::{self, record::GenotypeAllele::Unphased, Read};
use std::collections::{BTreeMap, HashMap};
use std::{path::PathBuf, str};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    hdf5_reader: hdf5::File,
    haplotype_variants: bcf::Reader,
    haplotype_calls: bcf::Reader,
    min_norm_counts: f64,
    outcsv: Option<PathBuf>,
}

impl Caller {
    pub fn call(&mut self) -> Result<()> {
        // Step 1: obtain kallisto estimates.
        let kallisto_estimates = KallistoEstimates::new(&self.hdf5_reader, self.min_norm_counts)?;
        dbg!(&kallisto_estimates);
        let haplotypes: Vec<String> = kallisto_estimates.keys().map(|x| x.to_string()).collect();
        let haplotype_variants = HaplotypeVariants::new(&mut self.haplotype_variants, &haplotypes)?;
        let haplotype_calls = HaplotypeCalls::new(&mut self.haplotype_calls)?;
        // Step 2: setup model.
        let model = Model::new(Likelihood::new(), Prior::new(), Posterior::new());

        //let universe = HaplotypeFractions::likely(&kallisto_estimates);
        let data = Data::new(
            kallisto_estimates.values().cloned().collect(),
            haplotype_variants,
            haplotype_calls,
        );

        // Step 3: calculate posteriors.
        //let m = model.compute(universe, &data);
        let m = model.compute_from_marginal(&Marginal::new(haplotypes.len()), &data);

        // Step 4: print TSV table with results
        // TODO use csv crate
        // Columns: posterior_prob, haplotype_a, haplotype_b, haplotype_c, ...
        // with each column after the first showing the fraction of the respective haplotype
        let mut posterior = m.event_posteriors();
        let mut wtr = csv::Writer::from_path(self.outcsv.as_ref().unwrap())?;
        let mut headers: Vec<_> = vec!["density".to_string(), "odds".to_string()];
        headers.extend(haplotypes);
        wtr.write_record(&headers)?;

        //write best record on top
        let mut records = Vec::new();
        let (haplotype_frequencies, best_density) = posterior.next().unwrap();
        let best_odds = 1;
        if best_density.exp() < 0.1 {
            records.push(format!("{:+.2e}", best_density.exp()));
        } else {
            records.push(best_density.exp().to_string());
        }
        records.push(best_odds.to_string());

        for frequency in haplotype_frequencies.iter() {
            if frequency < &NotNan::new(0.1).unwrap() {
                records.push(format!("{:+.2e}", NotNan::into_inner(*frequency)));
            } else {
                records.push(frequency.to_string())
            }
        }
        wtr.write_record(records)?;

        //write the rest of the records
        for (haplotype_frequencies, density) in posterior {
            let mut records = Vec::new();
            let odds = (density - best_density).exp();
            if density.exp() < 0.1 {
                records.push(format!("{:+.2e}", density.exp()));
            } else {
                records.push(density.exp().to_string());
            }

            if odds < 0.1 {
                records.push(format!("{:+.2e}", odds));
            } else {
                records.push(odds.to_string());
            }
            for frequency in haplotype_frequencies.iter() {
                if frequency < &NotNan::new(0.1).unwrap() {
                    records.push(format!("{:+.2e}", NotNan::into_inner(*frequency)));
                } else {
                    records.push(frequency.to_string())
                }
            }
            wtr.write_record(records)?;
        }
        Ok(())
    }
}

#[derive(Derefable, Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct Haplotype(#[deref] String);

#[derive(Debug, Clone)]
pub(crate) struct KallistoEstimate {
    pub count: f64,
    pub dispersion: f64,
}

#[derive(Debug, Clone, Derefable)]
pub(crate) struct KallistoEstimates(#[deref] HashMap<Haplotype, KallistoEstimate>);

impl KallistoEstimates {
    /// Generate new instance.
    pub(crate) fn new(hdf5_reader: &hdf5::File, min_norm_counts: f64) -> Result<Self> {
        let seqnames = Self::filter_seqnames(hdf5_reader, min_norm_counts)?;
        let ids = hdf5_reader
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::FixedAscii<255>>()?;
        let num_bootstraps = hdf5_reader.dataset("aux/num_bootstrap")?.read_1d::<i32>()?;
        let seq_length = hdf5_reader.dataset("aux/lengths")?.read_1d::<f64>()?;

        let mut estimates = HashMap::new();

        for seqname in seqnames {
            let index = ids.iter().position(|x| x.as_str() == seqname).unwrap();
            let mut bootstraps = Vec::new();
            for i in 0..num_bootstraps[0] {
                let dataset = hdf5_reader.dataset(&format!("bootstrap/bs{i}", i = i))?;
                let est_counts = dataset.read_1d::<f64>()?;
                let norm_counts = est_counts / &seq_length;
                let norm_counts = norm_counts[index];
                bootstraps.push(norm_counts);
            }
            //mean
            let sum = bootstraps.iter().sum::<f64>();
            let count = bootstraps.len();
            let m = sum / count as f64;

            //std dev
            let variance = bootstraps
                .iter()
                .map(|value| {
                    let diff = m - (*value as f64);
                    diff * diff
                })
                .sum::<f64>()
                / count as f64;
            let std = variance.sqrt();
            let t = std / m;

            //retrieval of mle
            let mle_dataset = hdf5_reader.dataset("est_counts")?.read_1d::<f64>()?;
            let mle_norm = mle_dataset / &seq_length; //normalized mle counts by length
            let m = mle_norm[index];

            estimates.insert(
                Haplotype(seqname.clone()),
                KallistoEstimate {
                    dispersion: t,
                    count: m,
                },
            );
        }

        Ok(KallistoEstimates(estimates))
    }

    //Return a vector of filtered seqnames according to --min-norm-counts.
    fn filter_seqnames(hdf5_reader: &hdf5::File, min_norm_counts: f64) -> Result<Vec<String>> {
        let est_counts = hdf5_reader.dataset("est_counts")?.read_1d::<f64>()?;
        let seq_length = hdf5_reader.dataset("aux/lengths")?.read_1d::<f64>()?; //these two variables arrays have the same length.
        let norm_counts = est_counts / seq_length;
        let mut indices = Vec::new();
        for (i, num) in norm_counts.iter().enumerate() {
            if num > &min_norm_counts {
                indices.push(i);
            }
        }
        let ids = hdf5_reader
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::FixedAscii<255>>()?;

        let mut filtered: Vec<String> = Vec::new();
        for i in indices {
            filtered.push(ids[i].to_string());
        }
        Ok(filtered)
    }
}

#[derive(Derefable, Debug, Clone, PartialEq, Eq, Hash, Ord, PartialOrd, new)]
pub(crate) struct VariantID(#[deref] i32);

#[derive(Derefable, Debug, Clone, PartialEq, Eq, PartialOrd)]
pub(crate) struct HaplotypeVariants(#[deref] BTreeMap<VariantID, (BitVec, BitVec)>);

impl HaplotypeVariants {
    pub(crate) fn new(
        haplotype_variants: &mut bcf::Reader,
        haplotypes: &Vec<String>,
    ) -> Result<Self> {
        let mut variant_records = BTreeMap::new();
        for record_result in haplotype_variants.records() {
            let record = record_result?;
            let variant_id: i32 = String::from_utf8(record.id())?.parse().unwrap();
            let header = record.header();
            let gts = record.genotypes()?; //genotypes of all samples
            let mut haplotype_indices = Vec::new();
            haplotypes.iter().for_each(|sample| {
                let position = header
                    .samples()
                    .iter()
                    .map(|x| str::from_utf8(x).unwrap())
                    .position(|x| x == sample)
                    .unwrap();
                haplotype_indices.push(position)
            });
            let mut haplotype_variants: BitVec<usize> = BitVec::new();
            haplotype_indices.iter().for_each(|sample_index| {
                for gta in gts.get(*sample_index).iter() {
                    haplotype_variants.push(*gta == Unphased(1));
                }
            });
            let mut haplotype_loci: BitVec<usize> = BitVec::new();
            let loci = record
                .format(b"C")
                .integer()
                .expect("Couldn't retrieve C field");
            haplotype_indices.iter().for_each(|sample_index| {
                if loci[*sample_index] == &[1] {}
                haplotype_loci.push(loci[*sample_index] == &[1]);
            });
            variant_records.insert(VariantID(variant_id), (haplotype_variants, haplotype_loci));
        }
        Ok(HaplotypeVariants(variant_records))
    }
}
#[derive(Debug, Clone, Derefable)]
pub(crate) struct AlleleFreqDist(#[deref] BTreeMap<AlleleFreq, f64>);

impl AlleleFreqDist {
    pub(crate) fn vaf_query(&self, vaf: AlleleFreq) -> LogProb {
        if self.contains_key(&vaf) {
            LogProb::from(PHREDProb(*self.get(&vaf).unwrap()))
        } else {
            let (x_0, y_0) = self.range(..vaf).next_back().unwrap();
            let (x_1, y_1) = self.range(vaf..).next().unwrap();
            let density = NotNan::new(*y_0).unwrap() + (vaf - *x_0) * (*y_1 - *y_0) / (*x_1 - *x_0); //calculation of density for given vaf by linear interpolation
            LogProb::from(PHREDProb(NotNan::into_inner(density)))
        }
    }
}

#[derive(Derefable, Debug, Clone)]
pub(crate) struct HaplotypeCalls(#[deref] BTreeMap<VariantID, AlleleFreqDist>);

impl HaplotypeCalls {
    pub(crate) fn new(haplotype_calls: &mut bcf::Reader) -> Result<Self> {
        let mut calls = BTreeMap::new();
        for record_result in haplotype_calls.records() {
            let record = record_result?;
            let afd_utf = record.format(b"AFD").string()?;
            if afd_utf[0] != b".".to_vec() {
                let variant_id: i32 = String::from_utf8(record.id())?.parse().unwrap();
                let afd = std::str::from_utf8(afd_utf[0]).unwrap();
                let mut vaf_density = BTreeMap::new();
                for pair in afd.split(',') {
                    let (vaf, density) = pair.split_once("=").unwrap();
                    let (vaf, density): (AlleleFreq, f64) =
                        (vaf.parse().unwrap(), density.parse().unwrap());
                    vaf_density.insert(vaf, density);
                }
                calls.insert(VariantID(variant_id), AlleleFreqDist(vaf_density));
            }
        }
        Ok(HaplotypeCalls(calls))
    }
}
// fn output_scientific(number: f64) -> {
//     if number < 0.1 {
//         records.push(format!("{:+.2e}", number.exp()));
//     } else {
//         records.push(number.exp().to_string());
//     }
// }

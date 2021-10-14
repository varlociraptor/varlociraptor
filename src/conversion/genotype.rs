use anyhow::Result;
use rust_htslib::bcf;
use rust_htslib::bcf::record::{GenotypeAllele, Numeric};
use rust_htslib::bcf::Read;

/// Decode PHRED scaled values to probabilities.
pub(crate) fn genotype() -> Result<()> {
    let mut inbcf = bcf::Reader::from_stdin()?;

    // setup output file
    let mut header = bcf::Header::from_template(inbcf.header());
    header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

    let mut outbcf = bcf::Writer::from_stdout(&header, false, bcf::Format::Bcf)?;

    for record in inbcf.records() {
        let mut record = record?;
        let vafs = record.format(b"AF").float()?;
        let dps = record.format(b"DP").integer()?;

        outbcf.translate(&mut record);
        let genotypes: Vec<_> = vafs
            .iter()
            .zip(dps.iter())
            .map(|(vaf, dp)| {
                let vaf = vaf[0];
                let dp = dp[0];
                if relative_eq!(vaf, 0.5) {
                    [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]
                } else if relative_eq!(vaf, 1.0) {
                    [GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)]
                } else if relative_eq!(vaf, 0.0) {
                    [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)]
                } else if !dp.is_missing() && dp > 0 {
                    // VAF is < 1.0 but not exactly 0.5. Still infer heterozygous genotype because
                    // it is the most likely case (in a subclone of the cells).
                    [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]
                } else {
                    // no observations
                    [
                        GenotypeAllele::UnphasedMissing,
                        GenotypeAllele::UnphasedMissing,
                    ]
                }
            })
            .flatten()
            .collect();
        record.push_genotypes(&genotypes)?;

        outbcf.write(&record)?;
    }

    Ok(())
}

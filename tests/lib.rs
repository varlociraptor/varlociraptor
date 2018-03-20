extern crate libprosic;
extern crate rust_htslib;
extern crate bio;
extern crate fern;
extern crate log;
extern crate itertools;
extern crate hyper;
extern crate flate2;

use std::fs;
use std::path::Path;
use std::str;
use std::{thread, time};
use std::io;
use std::process::Command;

use itertools::Itertools;
use rust_htslib::{bam,bcf};
use bio::stats::Prob;

use libprosic::model::{AlleleFreq, ContinuousAlleleFreqs, DiscreteAlleleFreqs};
use libprosic::constants;


fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}


fn cleanup_file(f: &str) {
    if Path::new(f).exists() {
        fs::remove_file(f).unwrap();
    }
}


pub fn setup_logger(test: &str) {
    let basedir = basedir(test);
    let logfile = format!("{}/debug.log", basedir);
    cleanup_file(&logfile);

    fern::Dispatch::new().level(log::LogLevelFilter::Debug)
                         .chain(fern::log_file(&logfile).unwrap())
                         .apply()
                         .unwrap();
    println!("Debug output can be found here: {}", logfile);
}


fn download_reference(chrom: &str) -> String {
    let reference = format!("tests/resources/{}.fa", chrom);
    if !Path::new(&reference).exists() {
        let client = hyper::Client::new();
        let res = client.get(
            &format!("http://hgdownload.cse.ucsc.edu/goldenpath/hg18/chromosomes/{}.fa.gz", chrom)
        ).send().unwrap();
        let mut reference_stream = flate2::read::GzDecoder::new(res).unwrap();
        let mut reference_file = fs::File::create(&reference).unwrap();

        io::copy(&mut reference_stream, &mut reference_file).unwrap();
    }
    assert!(Path::new(&reference).exists());
    if !Path::new(&(reference.clone() + ".fai")).exists() {
        Command::new("samtools").args(&["faidx", &reference])
                                .status()
                                .expect("failed to create fasta index");
    }

    reference
}


fn call_tumor_normal(test: &str, exclusive_end: bool, chrom: &str) {
    let reference = download_reference(chrom);

    //setup_logger(test);

    let basedir = basedir(test);

    let tumor_bam = format!("{}/tumor.bam", basedir);
    let normal_bam = format!("{}/normal.bam", basedir);

    let tumor_bam = bam::IndexedReader::from_path(&tumor_bam).unwrap();
    let normal_bam = bam::IndexedReader::from_path(&normal_bam).unwrap();

    let candidates = format!("{}/candidates.vcf", basedir);

    let output = format!("{}/calls.bcf", basedir);
    let observations = format!("{}/observations.tsv", basedir);
    cleanup_file(&output);
    cleanup_file(&observations);

    let insert_size = libprosic::InsertSize {
        mean: 312.0,
        sd: 15.0
    };
    let purity = 0.75;

    let tumor = libprosic::Sample::new(
        tumor_bam,
        2500,
        true,
        false,
        true,
        false,
        insert_size,
        libprosic::likelihood::LatentVariableModel::new(purity),
        constants::PROB_ILLUMINA_INS,
        constants::PROB_ILLUMINA_DEL,
        Prob(0.0),
        Prob(0.0),
        100
    );

    let normal = libprosic::Sample::new(
        normal_bam,
        2500,
        true,
        false,
        true,
        false,
        insert_size,
        libprosic::likelihood::LatentVariableModel::new(1.0),
        constants::PROB_ILLUMINA_INS,
        constants::PROB_ILLUMINA_DEL,
        Prob(0.0),
        Prob(0.0),
        100
    );

    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "germline".to_owned(),
            af_case: ContinuousAlleleFreqs::left_exclusive( 0.0..1.0 ),
            af_control: DiscreteAlleleFreqs::feasible(2).not_absent()
        },
        libprosic::call::pairwise::PairEvent {
            name: "somatic".to_owned(),
            af_case: ContinuousAlleleFreqs::left_exclusive( 0.0..1.0 ),
            af_control: DiscreteAlleleFreqs::absent()
        },
        libprosic::call::pairwise::PairEvent {
            name: "absent".to_owned(),
            af_case: ContinuousAlleleFreqs::left_exclusive( 0.0..0.0 ),
            af_control: DiscreteAlleleFreqs::absent()
        }
    ];

    let prior_model = libprosic::priors::FlatTumorNormalModel::new(2);

    let mut caller = libprosic::model::PairCaller::new(
        tumor,
        normal,
        prior_model
    );

    libprosic::call::pairwise::call::<
            _, _, _,
            libprosic::model::PairCaller<
                libprosic::model::ContinuousAlleleFreqs,
                libprosic::model::DiscreteAlleleFreqs,
                libprosic::model::priors::FlatTumorNormalModel
            >, _, _, _, _>
        (
            Some(&candidates),
            Some(&output),
            &reference,
            &events,
            &mut caller,
            false,
            false,
            Some(10000),
            Some(&observations),
            exclusive_end
        ).unwrap();
}

fn call_single_cell_bulk(test: &str, exclusive_end: bool, chrom: &str) {
    let reference = download_reference(chrom);

    //setup_logger(test);

    let basedir = basedir(test);

    let sc_bam_file = format!("{}/single_cell.bam", basedir);
    let bulk_bam_file = format!("{}/bulk.bam", basedir);

    let sc_bam = bam::IndexedReader::from_path(&sc_bam_file).unwrap();
    let bulk_bam = bam::IndexedReader::from_path(&bulk_bam_file).unwrap();

    let candidates = format!("{}/candidates.vcf", basedir);

    let output = format!("{}/calls.bcf", basedir);
    let observations = format!("{}/observations.tsv", basedir);
    cleanup_file(&output);
    cleanup_file(&observations);

    let insert_size = libprosic::InsertSize {
        mean: 312.0,
        sd: 15.0
    };

    let sc = libprosic::Sample::new(
        sc_bam,
        2500,
        true,
        true,
        true,
        false,
        insert_size,
        libprosic::likelihood::LatentVariableModel::with_single_sample(),
        constants::PROB_ILLUMINA_INS,
        constants::PROB_ILLUMINA_DEL,
        Prob(0.0),
        Prob(0.0),
        100
    );

    let bulk = libprosic::Sample::new(
        bulk_bam,
        2500,
        true,
        true,
        true,
        false,
        insert_size,
        libprosic::likelihood::LatentVariableModel::with_single_sample(),
        constants::PROB_ILLUMINA_INS,
        constants::PROB_ILLUMINA_DEL,
        Prob(0.0),
        Prob(0.0),
        100
    );

    // setup events: case = single cell; control = bulk
    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "hom_ref".to_owned(),
            af_case: DiscreteAlleleFreqs::absent(),
            af_control: ContinuousAlleleFreqs::right_exclusive( 0.0..0.5 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "ADO_to_ref".to_owned(),
            af_case: DiscreteAlleleFreqs::absent(),
            af_control: ContinuousAlleleFreqs::right_exclusive( 0.5..1.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "ADO_to_alt".to_owned(),
            af_case: DiscreteAlleleFreqs::new(vec![AlleleFreq(1.0)]),
            af_control: ContinuousAlleleFreqs::left_exclusive( 0.0..0.5 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "hom_alt".to_owned(),
            af_case: DiscreteAlleleFreqs::new(vec![AlleleFreq(1.0)]),
            af_control: ContinuousAlleleFreqs::left_exclusive( 0.5..1.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "err_alt".to_owned(),
            af_case: DiscreteAlleleFreqs::feasible(2).not_absent(),
            af_control: ContinuousAlleleFreqs::inclusive( 0.0..0.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "het".to_owned(),
            af_case: DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5)]),
            af_control: ContinuousAlleleFreqs::exclusive( 0.0..1.0 )
        },
        libprosic::call::pairwise::PairEvent {
            name: "err_ref".to_owned(),
            af_case: DiscreteAlleleFreqs::feasible(2).not_absent(),
            af_control: ContinuousAlleleFreqs::inclusive( 1.0..1.0 )
        }
    ];

    let prior_model = libprosic::priors::SingleCellBulkModel::new(2, 8, 100);

    let mut caller = libprosic::model::PairCaller::new(
        sc,
        bulk,
        prior_model
    );

    libprosic::call::pairwise::call::<
            _, _, _,
            libprosic::model::PairCaller<
                libprosic::model::DiscreteAlleleFreqs,
                libprosic::model::ContinuousAlleleFreqs,
                libprosic::model::priors::SingleCellBulkModel
            >, _, _, _, _>
        (
            Some(&candidates),
            Some(&output),
            &reference,
            &events,
            &mut caller,
            false,
            false,
            Some(10000),
            Some(&observations),
            exclusive_end
        ).unwrap();

    // sleep a second in order to wait for filesystem flushing
    thread::sleep(time::Duration::from_secs(1));
}

fn load_call(test: &str) -> bcf::Record {
    let basedir = basedir(test);

    let mut reader = bcf::Reader::from_path(format!("{}/calls.bcf", basedir)).unwrap();

    let mut calls = reader.records().map(|r| r.unwrap()).collect_vec();
    assert_eq!(calls.len(), 1, "unexpected number of calls");
    calls.pop().unwrap()
}


fn check_info_float(rec: &mut bcf::Record, tag: &[u8], truth: f32, maxerr: f32) {
    let p = rec.info(tag).float().unwrap().unwrap()[0];
    let err = (p - truth).abs();
    assert!(
        err <= maxerr,
        "{} error too high: value={}, truth={}, error={}",
        str::from_utf8(tag).unwrap(),
        p, truth, maxerr
    );
}


fn control_fdr_ev(test: &str) {
    let basedir = basedir(test);
    let mut calls = bcf::Reader::from_path(format!("{}/calls.matched.bcf", basedir)).unwrap();
    let output = format!("{}/thresholds.tsv", basedir);
    cleanup_file(&output);
    let mut writer = fs::File::create(&output).unwrap();
    libprosic::estimation::fdr::ev::control_fdr(
        &mut calls,
        &mut writer,
        &[ libprosic::SimpleEvent { name: "SOMATIC".to_owned() } ],
        &libprosic::model::VariantType::Deletion(Some(1..30))
    ).unwrap();
}


/// Test a Pindel call in a repeat region. It is a false positive that should be called as absent.
///
/// # Tweaks to make this work:
///
/// * we have to consider the full read for the PairHMM. This helps a lot in repeat regions.
#[test]
fn test1() {
    call_tumor_normal("test1", false, "chr1");
    let mut call = load_call("test1");

    check_info_float(&mut call, b"CASE_AF", 0.0, 0.15);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_ABSENT", 4.37, 0.05);
    // TODO weak enough for now, see if this can be improved
    check_info_float(&mut call, b"PROB_SOMATIC", 1.97, 0.05);
}


/// Test a Pindel call that is a somatic call in reality (case af: 0.125).
#[test]
fn test2() {
    call_tumor_normal("test2", false, "chr1");
    let mut call = load_call("test2");

    check_info_float(&mut call, b"CASE_AF", 0.125, 0.05);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_SOMATIC", 0.38, 0.01);
}


/// Test a Pindel call that is a germline call in reality (case af: 0.5, control af: 0.5).
#[test]
fn test3() {
    call_tumor_normal("test3", false, "chr1");
    let mut call = load_call("test3");

    check_info_float(&mut call, b"CASE_AF", 0.5, 0.05);
    check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.0);
    check_info_float(&mut call, b"PROB_GERMLINE", 7.0e-5, 7.0e-5);
}


/// Test a Pindel call (insertion) that is a somatic call in reality (case af: 0.042, control af: 0.0).
#[test]
fn test4() {
    call_tumor_normal("test4", false, "chr1");
    let mut call = load_call("test4");

    check_info_float(&mut call, b"CASE_AF", 0.042, 0.06);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_SOMATIC", 0.8, 0.5);
}


/// Test a Delly call in a repeat region. This should not be a somatic call.
#[test]
fn test5() {
    call_tumor_normal("test5", true, "chr1");
    let mut call = load_call("test5");
    check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.0);
    check_info_float(&mut call, b"PROB_SOMATIC", 41.4, 0.1);
}


/// Test a large deletion that should not be a somatic call. It seems to be absent, and has very
/// weak evidence.
#[test]
fn test6() {
    call_tumor_normal("test6", false, "chr16");
    let mut call = load_call("test6");
    check_info_float(&mut call, b"PROB_SOMATIC", 7.95, 0.01);
    //check_info_float(&mut call, b"PROB_ABSENT", 0.76, 0.01);
    //check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    //check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
}


/// Test a small Lancet deletion. It is a somatic call (AF=0.125) in reality.
#[test]
fn test7() {
    call_tumor_normal("test7", false, "chr1");
    let mut call = load_call("test7");
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"CASE_AF", 0.125, 0.12);
    check_info_float(&mut call, b"PROB_SOMATIC", 0.71, 0.01);
    check_info_float(&mut call, b"PROB_GERMLINE", 9.88, 0.01);
}


/// Test a Delly deletion. It is a germline call in reality.
#[test]
fn test8() {
    call_tumor_normal("test8", true, "chr2");
    let mut call = load_call("test8");
    check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.0);
    check_info_float(&mut call, b"PROB_GERMLINE", 0.0, 0.8);
}

/// Test a Delly deletion. It should not be a somatic call.
#[test]
fn test9() {
    call_tumor_normal("test9", true, "chr2");
    let mut call = load_call("test9");
    //check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.0);
    check_info_float(&mut call, b"PROB_SOMATIC", 2.98, 0.1);
}


/// Test a Lancet insertion. It seems to be a germline variant from venters genome. Evidence is
/// weak, but it should definitely not be called as somatic.
#[test]
fn test10() {
    call_tumor_normal("test10", false, "chr20");
    let mut call = load_call("test10");
    check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.0);
    check_info_float(&mut call, b"PROB_SOMATIC", 75.24, 0.05);
}


// A delly deletion that has very low coverage and very weak evidence. We cannot really infer
// something. However, this test is in here to ensure that such corner cases (a lot of -inf), do
// not cause panics.
#[test]
fn test11() {
    call_tumor_normal("test11", true, "chr2");
    load_call("test11");
}


/// A large lancet insertion that is not somatic, but likely a homozygous germline variant.
#[test]
fn test12() {
    call_tumor_normal("test12", false, "chr10");
    let mut call = load_call("test12");
    check_info_float(&mut call, b"CONTROL_AF", 1.0, 0.0);
    check_info_float(&mut call, b"CASE_AF", 1.0, 0.05);
}

/// A delly deletion that is a somatic mutation in reality (AF=0.33).
#[test]
fn test13() {
    call_tumor_normal("test13", true, "chr1");
    let mut call = load_call("test13");
    check_info_float(&mut call, b"PROB_SOMATIC", 0.00007, 0.00001);
    check_info_float(&mut call, b"CASE_AF", 0.33, 0.1);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
}

/// A delly deletion that is not somatic but a germline. However, there is a large bias
/// towards the ref allele in the normal sample. A reasonable explanation is a repeat structure
/// or amplification that projects reference allele reads on the variant location.
/// There is currently no way to avoid this, but an amplification factor could become another
/// latent variable in the model.
#[test]
#[ignore]
fn test14() {
    call_tumor_normal("test14", true, "chr15");
    let mut call = load_call("test14");
    check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.0);
}

/// A small lancet deletion that is a true and strong somatic variant (AF=1.0).
#[test]
fn test15() {
    call_tumor_normal("test15", false, "chr1");
    let mut call = load_call("test15");
    check_info_float(&mut call, b"PROB_SOMATIC", 1.2e-5, 1e-3);
    check_info_float(&mut call, b"CASE_AF", 1.0, 0.06);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
}

/// A large lancet deletion that is a true and strong somatic variant (AF=0.333).
#[test]
fn test16() {
    call_tumor_normal("test16", false, "chr1");
    let mut call = load_call("test16");
    check_info_float(&mut call, b"CASE_AF", 0.62, 0.3); // TODO this has a large bias
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_SOMATIC", 6.6e-5, 1e-3);
}

/// A delly call that is a false positive. It should be called as absent.
#[test]
fn test17() {
    call_tumor_normal("test17", true, "chr11");
    let mut call = load_call("test17");
    check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_ABSENT", 0.38, 0.01);
}

/// A large lancet deletion that is not somatic and a likely homozygous germline variant.
#[test]
fn test18() {
    call_tumor_normal("test18", false, "chr12");
    let mut call = load_call("test18");
    check_info_float(&mut call, b"CASE_AF", 1.0, 0.0);
    check_info_float(&mut call, b"CONTROL_AF", 1.0, 0.0);
    check_info_float(&mut call, b"PROB_GERMLINE", 0.0, 0.005);
}

/// A delly deletion that is not somatic but a heterozygous germline variant.
/// This needs handling of supplementary alignments when determining whether a fragment
/// is enclosing the variant.
#[test]
#[ignore]
fn test19() {
    call_tumor_normal("test19", true, "chr8");
    let mut call = load_call("test19");
    // TODO fix this
    check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.0);
}

/// A delly deletion that is not a somatic variant. It is in a highly repetetive region, and we
/// seem to be just unlucky that there is nothing visible in the normal sample.
#[test]
#[ignore]
fn test20() {
    call_tumor_normal("test20", true, "chr4");
    let mut call = load_call("test20");
    // TODO fix this
    check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.0);
    check_info_float(&mut call, b"PROB_GERMLINE", 0.19, 0.005);
}

/// A lancet insertion that is at the same place as a real somatic insertion, however
/// lancet calls a too short sequence. Ideally, we would call this as absent, but reads
/// are too small to do this properly.
#[test]
#[ignore]
fn test21() {
    call_tumor_normal("test21", false, "chr7");
    load_call("test21");
    assert!(false);
}

/// A manta deletion that is not a somatic variant. The deletion is also present in the normal,
/// but the allele frequency is too far away from 50% to be called as heterozygous.
#[test]
fn test22() {
    call_tumor_normal("test22", false, "chr18");
    let mut call = load_call("test22");
    check_info_float(&mut call, b"PROB_SOMATIC", 0.6, 0.01); // weak enough
    //check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.0);
}

/// Test a manta deletion that is not somatic. Seems to be absent.
#[test]
fn test23() {
    call_tumor_normal("test23", false, "chr14");
    let mut call = load_call("test23");
    check_info_float(&mut call, b"PROB_ABSENT", 3.65, 0.01);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
}

/// Test a small strelka deletion. It is not somatic, and should be called as absent.
#[test]
fn test24() {
    call_tumor_normal("test24", false, "chr6");
    let mut call = load_call("test24");
    check_info_float(&mut call, b"PROB_SOMATIC", 5.3, 0.01);
    check_info_float(&mut call, b"PROB_ABSENT", 2.13, 0.01);
    check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
}


#[test]
fn test_fdr_ev1() {
    control_fdr_ev("test_fdr_ev_1");
    // TODO add a reasonable assertion
}


#[test]
fn test_fdr_ev2() {
    control_fdr_ev("test_fdr_ev_2");
    // TODO add a reasonable assertion
}

#[test]
fn test_sc_bulk() {
    call_single_cell_bulk("test_sc_bulk", true, "chr1");
    let mut call = load_call("test_sc_bulk");
    check_info_float(&mut call, b"CONTROL_AF", 0.0285714, 0.0000001);
    check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_HET", 12.5142, 0.0001);
}

#[test]
fn test_sc_bulk_hom_ref() {
    call_single_cell_bulk("test_sc_bulk_hom_ref", true, "chr1");
    let mut call = load_call("test_sc_bulk_hom_ref");
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_HOM_REF", 0.219712, 0.000001);
}

#[test]
fn test_sc_bulk_indel() {
    call_single_cell_bulk("test_sc_bulk_indel", true, "chr1");
    let mut call = load_call("test_sc_bulk_indel");
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_HET", 11.26, 0.01);
}

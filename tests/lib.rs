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

use itertools::Itertools;
use rust_htslib::{bam,bcf};
use bio::stats::Prob;

use libprosic::model::AlleleFreq;
use libprosic::constants;


fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}


fn cleanup_file(f: &str) {
    if Path::new(f).exists() {
        fs::remove_file(f).unwrap();
    }
}


fn setup_logger(test: &str) {
    let basedir = basedir(test);
    let logfile = format!("{}/debug.log", basedir);
    cleanup_file(&logfile);

    fern::Dispatch::new().level(log::LogLevelFilter::Debug)
                         .chain(fern::log_file(&logfile).unwrap())
                         .apply()
                         .unwrap();
    println!("Debug output can be found here: {}", logfile);
}


fn download_reference() -> &'static str {
    let reference = "tests/resources/chr1.fa";
    if !Path::new(reference).exists() {
        let client = hyper::Client::new();
        let res = client.get(
            "http://hgdownload.cse.ucsc.edu/goldenpath/hg18/chromosomes/chr1.fa.gz"
        ).send().unwrap();
        let mut reference_stream = flate2::read::GzDecoder::new(res).unwrap();
        let mut reference_file = fs::File::create(reference).unwrap();

        io::copy(&mut reference_stream, &mut reference_file).unwrap();
    }
    reference
}


fn call_tumor_normal(test: &str, exclusive_end: bool) {
    let reference = download_reference();
    assert!(Path::new(reference).exists());
    assert!(Path::new(&(reference.to_owned() + ".fai")).exists());

    //setup_logger(test);

    let basedir = basedir(test);

    let tumor_bam = format!("{}/tumor.bam", basedir);
    let normal_bam = format!("{}/normal.bam", basedir);

    let tumor_bam = bam::IndexedReader::from_path(&tumor_bam).unwrap();
    let normal_bam = bam::IndexedReader::from_path(&normal_bam).unwrap();

    let candidates = format!("{}/candidates.vcf", basedir);
    let reference = "tests/resources/chr1.fa";

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
        true,
        true,
        false,
        insert_size,
        libprosic::likelihood::LatentVariableModel::new(purity),
        Prob(0.0),
        constants::PROB_ILLUMINA_INS,
        constants::PROB_ILLUMINA_DEL,
        Prob(0.0),
        Prob(0.0),
        25,
        25
    );

    let normal = libprosic::Sample::new(
        normal_bam,
        2500,
        true,
        true,
        true,
        false,
        insert_size,
        libprosic::likelihood::LatentVariableModel::new(1.0),
        Prob(0.0),
        constants::PROB_ILLUMINA_INS,
        constants::PROB_ILLUMINA_DEL,
        Prob(0.0),
        Prob(0.0),
        25,
        25
    );

    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "germline".to_owned(),
            af_case: AlleleFreq(0.0)..AlleleFreq(1.0),
            af_control: vec![AlleleFreq(0.5), AlleleFreq(1.0)]
        },
        libprosic::call::pairwise::PairEvent {
            name: "somatic".to_owned(),
            af_case: AlleleFreq(0.0)..AlleleFreq(1.0),
            af_control: vec![AlleleFreq(0.0)]
        },
        libprosic::call::pairwise::PairEvent {
            name: "absent".to_owned(),
            af_case: AlleleFreq(0.0)..AlleleFreq(0.0),
            af_control: vec![AlleleFreq(0.0)]
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

    // sleep a second in order to wait for filesystem flushing
    thread::sleep(time::Duration::from_secs(1));
}


fn load_call(test: &str) -> bcf::Record {
    let basedir = basedir(test);

    let reader = bcf::Reader::from_path(format!("{}/calls.bcf", basedir)).unwrap();

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
        &libprosic::SimpleEvent { name: "SOMATIC".to_owned() },
        &libprosic::model::VariantType::Deletion(Some(1..30))
    ).unwrap();
}


/// Test a Pindel call in a repeat region. It is a false positive that should be absent instead.
/// But so far it is reported as somatic (however, at least with a weak probability).
#[test]
// TODO remove ignore once fixed
#[ignore]
fn test1() {
    call_tumor_normal("test1", false);
    let mut call = load_call("test1");

    check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
}


/// Test a Pindel call that is a somatic call in reality (case af: 0.125).
#[test]
fn test2() {
    call_tumor_normal("test2", false);
    let mut call = load_call("test2");

    check_info_float(&mut call, b"CASE_AF", 0.125, 0.05);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_SOMATIC", 1.0e-7, 1.0e-7);
}


/// Test a Pindel call that is a germline call in reality (case af: 0.5, control af: 0.5).
#[test]
fn test3() {
    call_tumor_normal("test3", false);
    let mut call = load_call("test3");

    check_info_float(&mut call, b"CASE_AF", 0.5, 0.0);
    check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.0);
    check_info_float(&mut call, b"PROB_GERMLINE", 8.0e-29, 8.0e-29);
}


/// Test a Pindel call (insertion) that is a somatic call in reality (case af: 0.042, control af: 0.0).
#[test]
fn test4() {
    call_tumor_normal("test4", false);
    let mut call = load_call("test4");

    check_info_float(&mut call, b"CASE_AF", 0.042, 0.1);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_SOMATIC", 0.01, 0.06);
}


/// Test a Delly call in a repeat region. It is a false positive that should be absent instead.
/// But so far it is reported as somatic (however, at least with a weak probability).
#[test]
// TODO remove ignore once fixed
#[ignore]
fn test5() {
    call_tumor_normal("test5", true);
    let mut call = load_call("test5");
    check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
}


#[test]
fn test_fdr_ev() {
    control_fdr_ev("test_fdr_ev_1");
    // TODO add a reasonable assertion
}

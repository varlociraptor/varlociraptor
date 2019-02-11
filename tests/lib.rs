extern crate bio;
extern crate csv;
extern crate fern;
extern crate itertools;
extern crate libprosic;
extern crate log;
extern crate rust_htslib;
extern crate serde_json;

use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::str;
use std::{thread, time};

use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use rust_htslib::bcf::Read;
use rust_htslib::{bam, bcf};

use libprosic::constants;
use libprosic::model::{AlleleFreq, ContinuousAlleleFreqs, DiscreteAlleleFreqs};

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

    fern::Dispatch::new()
        .level(log::LogLevelFilter::Debug)
        .chain(fern::log_file(&logfile).unwrap())
        .apply()
        .unwrap();
    println!("Debug output can be found here: {}", logfile);
}

fn download_reference(chrom: &str, build: &str) -> PathBuf {
    let p = format!("tests/resources/{}/{}.fa", build, chrom);
    let reference = Path::new(&p);
    fs::create_dir_all(reference.parent().unwrap()).unwrap();
    if !reference.exists() {
        let url = if build.starts_with("hg") {
            format!(
                "http://hgdownload.cse.ucsc.edu/goldenpath/{}/chromosomes/{}.fa.gz",
                build, chrom
            )
        } else if build.starts_with("GRCh") {
            format!(
                // use soft-masked DNA sequence because we need it to omit variants within repeats.
                "ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.{}.dna_sm.chromosome.{}.fa.gz",
                build, chrom
            )
        } else {
            panic!("invalid genome build: {}", build);
        };

        let curl = Command::new("curl")
            .arg(&url)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn().unwrap();
        let mut gzip = Command::new("gzip")
            .arg("-d")
            .stdin(curl.stdout.unwrap())
            .stdout(Stdio::piped())
            .spawn().unwrap();
        let mut reference_file = fs::File::create(&reference).unwrap();
        io::copy(gzip.stdout.as_mut().unwrap(), &mut reference_file).unwrap();
    }

    assert!(Path::new(&reference).exists());
    if !reference.with_extension("fa.fai").exists() {
        Command::new("samtools")
            .args(&["faidx", reference.to_str().unwrap()])
            .status()
            .expect("failed to create fasta index");
    }

    reference.to_path_buf()
}

fn call_tumor_normal(test: &str, exclusive_end: bool, purity: f64, chrom: &str, build: &str) {
    let reference = download_reference(chrom, build);

    //setup_logger(test);

    let basedir = basedir(test);

    let tumor_bam_path = format!("{}/tumor.bam", basedir);
    let normal_bam_path = format!("{}/normal.bam", basedir);

    let tumor_bam = bam::IndexedReader::from_path(&tumor_bam_path).unwrap();
    let normal_bam = bam::IndexedReader::from_path(&normal_bam_path).unwrap();

    let candidates = format!("{}/candidates.vcf", basedir);
    let alignment_properties_def = format!("{}/alignment_properties.json", basedir);
    let alignment_properties_tumor_def = format!("{}/alignment_properties_tumor.json", basedir);
    let alignment_properties_normal_def = format!("{}/alignment_properties_normal.json", basedir);

    let output = format!("{}/calls.bcf", basedir);
    let observations = format!("{}/observations.tsv", basedir);
    cleanup_file(&output);
    cleanup_file(&observations);

    let (alignment_properties_tumor, alignment_properties_normal) = if Path::new(&alignment_properties_tumor_def).exists() {
        (
            serde_json::from_reader(fs::File::open(alignment_properties_tumor_def).unwrap()).unwrap(),
            serde_json::from_reader(fs::File::open(alignment_properties_normal_def).unwrap()).unwrap()
        )
    } else if Path::new(&alignment_properties_def).exists() {
        let a = serde_json::from_reader(fs::File::open(alignment_properties_def).unwrap()).unwrap();
        (a, a)
    } else {
        let mut bam = bam::Reader::from_path("tests/resources/tumor-first30000.bam").unwrap();
        //let mut bam = bam::Reader::from_path(&normal_bam_path).unwrap();
        let a = libprosic::AlignmentProperties::estimate(&mut bam).unwrap();
        (a, a)
    };
    println!("{:?}", alignment_properties_tumor);
    println!("{:?}", alignment_properties_normal);

    let tumor = libprosic::Sample::new(
        tumor_bam,
        true,
        alignment_properties_tumor,
        libprosic::likelihood::LatentVariableModel::new(purity),
        constants::PROB_ILLUMINA_INS,
        constants::PROB_ILLUMINA_DEL,
        Prob(0.0),
        Prob(0.0),
        100,
        500,
        &[],
    );

    let normal = libprosic::Sample::new(
        normal_bam,
        true,
        alignment_properties_normal,
        libprosic::likelihood::LatentVariableModel::new(1.0),
        constants::PROB_ILLUMINA_INS,
        constants::PROB_ILLUMINA_DEL,
        Prob(0.0),
        Prob(0.0),
        100,
        500,
        &[],
    );

    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "germline_het".to_owned(),
            af_case: ContinuousAlleleFreqs::inclusive(0.0..1.0),
            af_control: ContinuousAlleleFreqs::singleton(0.5),
        },
        libprosic::call::pairwise::PairEvent {
            name: "germline_hom".to_owned(),
            af_case: ContinuousAlleleFreqs::inclusive(0.0..1.0),
            af_control: ContinuousAlleleFreqs::singleton(1.0),
        },
        libprosic::call::pairwise::PairEvent {
            name: "somatic_tumor".to_owned(),
            af_case: ContinuousAlleleFreqs::left_exclusive(0.0..1.0).min_observations(2),
            af_control: ContinuousAlleleFreqs::absent(),
        },
        libprosic::call::pairwise::PairEvent {
            name: "somatic_normal".to_owned(),
            af_case: ContinuousAlleleFreqs::left_exclusive(0.0..1.0).min_observations(2),
            af_control: ContinuousAlleleFreqs::exclusive(0.0..0.5).min_observations(2),
        },
        libprosic::call::pairwise::PairEvent {
            name: "absent".to_owned(),
            af_case: ContinuousAlleleFreqs::absent(),
            af_control: ContinuousAlleleFreqs::absent(),
        },
    ];

    let prior_model = libprosic::priors::FlatTumorNormalModel::new(2);
    //let prior_model = libprosic::priors::TumorNormalModel::new(2, 3000.0, 0.01, 0.01, 3500000000, Prob(0.001));

    let mut caller = libprosic::model::PairCaller::new(tumor, normal, prior_model);

    libprosic::call::pairwise::call::<
        _,
        _,
        _,
        libprosic::model::PairCaller<
            libprosic::model::ContinuousAlleleFreqs,
            libprosic::model::ContinuousAlleleFreqs,
            libprosic::model::priors::FlatTumorNormalModel,
        >,
        _,
        _,
        _,
        _,
    >(
        Some(&candidates),
        Some(&output),
        &reference.to_str().unwrap(),
        &events,
        &mut caller,
        false,
        false,
        Some(10000),
        Some(&observations),
        exclusive_end,
    )
    .unwrap();
}

fn call_single_cell_bulk(test: &str, exclusive_end: bool, chrom: &str, build: &str) {
    let reference = download_reference(chrom, build);

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
        sd: 15.0,
    };
    let alignment_properties = libprosic::AlignmentProperties::default(insert_size);

    let sc = libprosic::Sample::new(
        sc_bam,
        true,
        alignment_properties,
        libprosic::likelihood::LatentVariableModel::with_single_sample(),
        constants::PROB_ILLUMINA_INS,
        constants::PROB_ILLUMINA_DEL,
        Prob(0.0),
        Prob(0.0),
        100,
        500,
        &[],
    );

    let bulk = libprosic::Sample::new(
        bulk_bam,
        true,
        alignment_properties,
        libprosic::likelihood::LatentVariableModel::with_single_sample(),
        constants::PROB_ILLUMINA_INS,
        constants::PROB_ILLUMINA_DEL,
        Prob(0.0),
        Prob(0.0),
        100,
        500,
        &[],
    );

    // setup events: case = single cell; control = bulk
    let events = [
        libprosic::call::pairwise::PairEvent {
            name: "hom_ref".to_owned(),
            af_case: DiscreteAlleleFreqs::absent(),
            af_control: ContinuousAlleleFreqs::right_exclusive(0.0..0.5),
        },
        libprosic::call::pairwise::PairEvent {
            name: "ADO_to_ref".to_owned(),
            af_case: DiscreteAlleleFreqs::absent(),
            af_control: ContinuousAlleleFreqs::right_exclusive(0.5..1.0),
        },
        libprosic::call::pairwise::PairEvent {
            name: "ADO_to_alt".to_owned(),
            af_case: DiscreteAlleleFreqs::new(vec![AlleleFreq(1.0)]),
            af_control: ContinuousAlleleFreqs::left_exclusive(0.0..0.5),
        },
        libprosic::call::pairwise::PairEvent {
            name: "hom_alt".to_owned(),
            af_case: DiscreteAlleleFreqs::new(vec![AlleleFreq(1.0)]),
            af_control: ContinuousAlleleFreqs::left_exclusive(0.5..1.0),
        },
        libprosic::call::pairwise::PairEvent {
            name: "err_alt".to_owned(),
            af_case: DiscreteAlleleFreqs::feasible(2).not_absent(),
            af_control: ContinuousAlleleFreqs::inclusive(0.0..0.0),
        },
        libprosic::call::pairwise::PairEvent {
            name: "het".to_owned(),
            af_case: DiscreteAlleleFreqs::new(vec![AlleleFreq(0.5)]),
            af_control: ContinuousAlleleFreqs::exclusive(0.0..1.0),
        },
        libprosic::call::pairwise::PairEvent {
            name: "err_ref".to_owned(),
            af_case: DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0), AlleleFreq(0.5)]),
            af_control: ContinuousAlleleFreqs::inclusive(1.0..1.0),
        },
    ];

    let prior_model = libprosic::priors::SingleCellBulkModel::new(2, 8, 100);

    let mut caller = libprosic::model::PairCaller::new(sc, bulk, prior_model);

    libprosic::call::pairwise::call::<
        _,
        _,
        _,
        libprosic::model::PairCaller<
            libprosic::model::DiscreteAlleleFreqs,
            libprosic::model::ContinuousAlleleFreqs,
            libprosic::model::priors::SingleCellBulkModel,
        >,
        _,
        _,
        _,
        _,
    >(
        Some(&candidates),
        Some(&output),
        &reference,
        &events,
        &mut caller,
        false,
        false,
        Some(10000),
        Some(&observations),
        exclusive_end,
    )
    .unwrap();

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
        p,
        truth,
        maxerr
    );
}

fn check_info_float_at_most(rec: &mut bcf::Record, tag: &[u8], threshold: f32) {
    let p = rec.info(tag).float().unwrap().unwrap()[0];
    assert!(
        p <= threshold,
        "{} too large: value={}, threshold={}",
        str::from_utf8(tag).unwrap(),
        p,
        threshold,
    );
}

fn check_info_float_at_least(rec: &mut bcf::Record, tag: &[u8], threshold: f32) {
    let p = rec.info(tag).float().unwrap().unwrap()[0];
    assert!(
        p >= threshold,
        "{} too small: value={}, threshold={}",
        str::from_utf8(tag).unwrap(),
        p,
        threshold,
    );
}

fn assert_call_number(test: &str, expected_calls: usize) {
    let basedir = basedir(test);

    let mut reader = bcf::Reader::from_path(format!("{}/calls.filtered.bcf", basedir)).unwrap();

    let calls = reader.records().map(|r| r.unwrap()).collect_vec();
    // allow one more or less, in order to be robust to numeric fluctuations
    assert!(
        (calls.len() as i32 - expected_calls as i32).abs() <= 1,
        "unexpected number of calls"
    );
}

fn control_fdr_ev(test: &str, event_str: &str, alpha: f64) {
    let basedir = basedir(test);
    let output = format!("{}/calls.filtered.bcf", basedir);
    cleanup_file(&output);
    libprosic::filtration::fdr::control_fdr(
        &format!("{}/calls.matched.bcf", basedir),
        Some(&output),
        &[libprosic::SimpleEvent {
            name: event_str.to_owned(),
        }],
        &libprosic::model::VariantType::Deletion(Some(1..30)),
        LogProb::from(Prob(alpha)),
    )
    .unwrap();
}

/// Test a Pindel call in a repeat region. It is either germline or absent, and could be called either
/// as deletion or insertion due to the special repeat structure here.
#[test]
fn test01() {
    call_tumor_normal("test1", false, 0.75, "chr1", "hg18");
    let mut call = load_call("test1");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 17.8);
}

/// Test a Pindel call that is a somatic call in reality (case af: 0.125).
#[test]
fn test02() {
    call_tumor_normal("test2", false, 0.75, "chr1", "hg18");
    let mut call = load_call("test2");

    check_info_float(&mut call, b"CASE_AF", 0.125, 0.07);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    // There are some fragments with increased insert size that let prosic think there is a bit of
    // evidence for having a somatic normal call. This makes the probability for somatic tumor weak.
    // TODO We could avoid this by raising the lower AF bound for somatic normal calls.
    check_info_float_at_most(&mut call, b"PROB_SOMATIC_TUMOR", 0.12);
}

/// Test a Pindel call that is a germline call in reality (case af: 0.5, control af: 0.5).
#[test]
fn test03() {
    call_tumor_normal("test3", false, 0.75, "chr1", "hg18");
    let mut call = load_call("test3");

    check_info_float(&mut call, b"CASE_AF", 0.5, 0.04);
    check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.051);
    check_info_float_at_most(&mut call, b"PROB_GERMLINE_HET", 0.8);
}

/// Test a Pindel call (insertion) that is a somatic call in reality (case af: 0.042, control af: 0.0).
#[test]
fn test04() {
    call_tumor_normal("test4", false, 0.75, "chr1", "hg18");
    let mut call = load_call("test4");

    check_info_float(&mut call, b"CASE_AF", 0.042, 0.1);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float_at_most(&mut call, b"PROB_SOMATIC_TUMOR", 0.06);
}

/// Test a Delly call in a repeat region. This should not be a somatic call.
#[test]
fn test05() {
    call_tumor_normal("test5", true, 0.75, "chr1", "hg18");
    let mut call = load_call("test5");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 435.0);
}

/// Test a large deletion that should not be a somatic call. It seems to be germline but a bit
/// unclear because of being in a repetetive region.
#[test]
fn test06() {
    call_tumor_normal("test6", false, 0.75, "chr16", "hg18");
    let mut call = load_call("test6");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 12.0);
}

/// Test a small Lancet deletion. It is a somatic call (AF=0.125) in reality.
#[test]
fn test07() {
    call_tumor_normal("test7", false, 0.75, "chr1", "hg18");
    let mut call = load_call("test7");
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"CASE_AF", 0.125, 0.06);
    check_info_float_at_most(&mut call, b"PROB_SOMATIC_TUMOR", 0.13);
}

/// Test a Delly deletion. It is a germline call in reality.
#[test]
fn test08() {
    call_tumor_normal("test8", true, 0.75, "chr2", "hg18");
    let mut call = load_call("test8");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 1232.0);
}

/// Test a Delly deletion. It should not be a somatic call.
#[test]
fn test09() {
    call_tumor_normal("test9", true, 0.75, "chr2", "hg18");
    let mut call = load_call("test9");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 110.0);
}

/// Test a Lancet insertion. It seems to be a germline variant from venters genome. Evidence is
/// weak, but it should definitely not be called as somatic.
#[test]
fn test10() {
    call_tumor_normal("test10", false, 0.75, "chr20", "hg18");
    let mut call = load_call("test10");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 1089.0);
}

// A delly deletion that has very low coverage and very weak evidence. We cannot really infer
// something. However, this test is in here to ensure that such corner cases (a lot of -inf), do
// not cause panics.
#[test]
fn test11() {
    call_tumor_normal("test11", true, 0.75, "chr2", "hg18");
    load_call("test11");
}

/// A large lancet insertion that is not somatic, but likely a homozygous germline variant.
#[test]
fn test12() {
    call_tumor_normal("test12", false, 0.75, "chr10", "hg18");
    let mut call = load_call("test12");
    check_info_float(&mut call, b"CONTROL_AF", 1.0, 0.0);
    check_info_float(&mut call, b"CASE_AF", 1.0, 0.0);
    check_info_float(&mut call, b"PROB_GERMLINE_HOM", 0.0, 0.01);
}

/// A delly deletion that is a somatic mutation in reality (AF=0.33).
#[test]
fn test13() {
    call_tumor_normal("test13", true, 0.75, "chr1", "hg18");
    let mut call = load_call("test13");
    check_info_float_at_most(&mut call, b"PROB_SOMATIC_TUMOR", 0.17);
    check_info_float(&mut call, b"CASE_AF", 0.33, 0.19);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
}

/// A delly deletion that is not somatic but a germline. However, there is a large bias
/// towards the ref allele in the normal sample. A reasonable explanation is a repeat structure
/// or amplification that projects reference allele reads on the variant location.
/// There is currently no way to avoid this, but an amplification factor could become another
/// latent variable in the model.
#[test]
fn test14() {
    call_tumor_normal("test14", true, 0.75, "chr15", "hg18");
    let mut call = load_call("test14");
    check_info_float(&mut call, b"CONTROL_AF", 0.5, 0.4);
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 349.0)
}

/// A small lancet deletion that is a true and strong somatic variant (AF=1.0).
#[test]
fn test15() {
    call_tumor_normal("test15", false, 0.75, "chr1", "hg18");
    let mut call = load_call("test15");
    check_info_float_at_most(&mut call, b"PROB_SOMATIC_TUMOR", 0.08);
    check_info_float(&mut call, b"CASE_AF", 1.0, 0.06);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
}

/// A large lancet deletion that is a true and strong somatic variant (AF=0.333).
#[test]
fn test16() {
    call_tumor_normal("test16", false, 0.75, "chr1", "hg18");
    let mut call = load_call("test16");
    check_info_float(&mut call, b"CASE_AF", 0.333, 0.12);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float_at_most(&mut call, b"PROB_SOMATIC_TUMOR", 0.12);
}

/// A delly call that is a false positive. It should be called as absent.
#[test]
fn test17() {
    call_tumor_normal("test17", true, 0.75, "chr11", "hg18");
    let mut call = load_call("test17");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 4.8);
    check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_ABSENT", 0.06, 0.01);
}

/// A large lancet deletion that is not somatic and a likely homozygous germline variant.
#[test]
fn test18() {
    call_tumor_normal("test18", false, 0.75, "chr12", "hg18");
    let mut call = load_call("test18");
    check_info_float(&mut call, b"CASE_AF", 1.0, 0.0);
    check_info_float(&mut call, b"CONTROL_AF", 1.0, 0.0);
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 2681.0);
}

/// A delly deletion that is not somatic but a heterozygous germline variant.
/// This needs handling of supplementary alignments when determining whether a fragment
/// is enclosing the variant.
#[test]
fn test19() {
    call_tumor_normal("test19", true, 0.75, "chr8", "hg18");
    let mut call = load_call("test19");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 446.0);
}

/// A delly deletion that is not a somatic variant but germline. It is in a highly repetetive
/// region, which causes a bias on MAPQ for variant fragments. However, it works when considering
/// MAPQ to be binary.
#[test]
fn test20() {
    call_tumor_normal("test20", true, 0.75, "chr4", "hg18");
    let mut call = load_call("test20");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 432.0);
}

/// A lancet insertion that is at the same place as a real somatic insertion, however
/// lancet calls a too short sequence. Ideally, we would call this as absent, but reads
/// are too small to do this properly.
#[test]
#[ignore]
fn test21() {
    call_tumor_normal("test21", false, 0.75, "chr7", "hg18");
    load_call("test21");
    assert!(false);
}

/// A manta deletion that is a germline variant. It seems to be homozygous when looking at IGV though.
#[test]
fn test22() {
    call_tumor_normal("test22", false, 0.75, "chr18", "hg18");
    let mut call = load_call("test22");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 152.0);
}

/// Test a manta deletion that is not somatic. It might be a repeat artifact or a heterozygous
/// germline variant.
#[test]
fn test23() {
    call_tumor_normal("test23", false, 0.75, "chr14", "hg18");
    let mut call = load_call("test23");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 6.6);
}

/// Test a small strelka deletion that is not somatic.
#[test]
fn test24() {
    call_tumor_normal("test24", false, 0.75, "chr6", "hg18");
    let mut call = load_call("test24");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 162.0);
}

/// Test a small lancet deletion that is a clear germline variant.
#[test]
fn test25() {
    call_tumor_normal("test25", false, 0.75, "chr11", "hg18");
    let mut call = load_call("test25");
    check_info_float(&mut call, b"CASE_AF", 1.0, 0.0);
    check_info_float(&mut call, b"CONTROL_AF", 1.0, 0.0);
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 2648.0);
}

/// Test a delly deletion (on real data) that is likely a repeat artifact.
#[test]
fn test26() {
    call_tumor_normal("test26", true, 1.0, "1", "GRCh38");
    let mut call = load_call("test26");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 4.6);
}

/// Test a delly deletion that is not a somatic variant. It is likely absent.
#[test]
fn test27() {
    call_tumor_normal("test27", true, 0.75, "chr10", "hg18");
    let mut call = load_call("test27");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 14.0);
}

/// Test a delly deletion that is not a somatic variant. It is likely absent.
#[test]
fn test28() {
    call_tumor_normal("test28", true, 0.75, "chr5", "hg18");
    let mut call = load_call("test28");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 10.6);
}

/// Test a delly deletion that is likely an artifact of a repeat region.
#[test]
fn test29() {
    call_tumor_normal("test29", true, 1.0, "1", "GRCh38");
    let mut call = load_call("test29");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 54.0);
    check_info_float_at_most(&mut call, b"PROB_ABSENT", 0.26);
}

/// Test a delly deletion that is likely a germline variant.
#[test]
fn test30() {
    call_tumor_normal("test30", true, 1.0, "2", "GRCh38");
    let mut call = load_call("test30");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 2296.0);
}

/// Test a delly deletion that is not a somatic variant.
#[test]
fn test31() {
    call_tumor_normal("test31", true, 1.0, "1", "GRCh38");
    let mut call = load_call("test31");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 11.08);
}

/// Test a delly deletion that is not a somatic variant.
#[test]
fn test32() {
    call_tumor_normal("test32", true, 1.0, "1", "GRCh38");
    let mut call = load_call("test32");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 1.1);
}

/// Test a delly deletion that is not a somatic variant.
#[test]
fn test33() {
    call_tumor_normal("test33", true, 1.0, "1", "GRCh38");
    let mut call = load_call("test33");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 3.13);
}

/// Test a delly deletion that is not a somatic variant.
#[test]
fn test34() {
    call_tumor_normal("test34", true, 1.0, "18", "GRCh38");
    let mut call = load_call("test34");
    check_info_float_at_least(&mut call, b"PROB_SOMATIC_TUMOR", 3.13);
}

#[test]
fn test_fdr_ev1() {
    control_fdr_ev("test_fdr_ev_1", "SOMATIC", 0.05);
    //assert_call_number("test_fdr_ev_1", 974);
}

#[test]
fn test_fdr_ev2() {
    control_fdr_ev("test_fdr_ev_2", "SOMATIC", 0.05);
    assert_call_number("test_fdr_ev_2", 950);
}

/// same test, but low alpha
#[test]
fn test_fdr_ev3() {
    control_fdr_ev("test_fdr_ev_3", "ABSENT", 0.001);
    assert_call_number("test_fdr_ev_3", 0);
}

#[test]
fn test_fdr_ev4() {
    control_fdr_ev("test_fdr_ev_4", "SOMATIC_TUMOR", 0.05);
    assert_call_number("test_fdr_ev_4", 0);
}

#[test]
fn test_sc_bulk() {
    call_single_cell_bulk("test_sc_bulk", true, "chr1", "hg18");
    let mut call = load_call("test_sc_bulk");
    check_info_float(&mut call, b"CONTROL_AF", 0.0285714, 0.0000001);
    check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_HET", 12.5142, 0.0001);
}

#[test]
fn test_sc_bulk_hom_ref() {
    call_single_cell_bulk("test_sc_bulk_hom_ref", true, "chr1", "hg18");
    let mut call = load_call("test_sc_bulk_hom_ref");
    check_info_float(&mut call, b"CONTROL_AF", 0.0, 0.0);
    check_info_float(&mut call, b"CASE_AF", 0.0, 0.0);
    check_info_float(&mut call, b"PROB_HOM_REF", 0.219712, 0.000001);
}

#[test]
fn test_sc_bulk_indel() {
    call_single_cell_bulk("test_sc_bulk_indel", true, "chr1", "hg18");
    let mut call = load_call("test_sc_bulk_indel");
    check_info_float(&mut call, b"CONTROL_AF", 0.13, 0.01);
    check_info_float(&mut call, b"CASE_AF", 0.5, 0.0);
    println!(
        "PROB_HOM_REF: {}",
        call.info(b"PROB_HOM_REF").float().unwrap().unwrap()[0]
    );
    check_info_float_at_least(&mut call, b"PROB_HOM_REF", 3.7);
}

use std::fs;
use std::path::Path;
use std::str;

use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use rust_htslib::bcf;
use rust_htslib::bcf::Read as BCFRead;

mod common;

use common::load_testcase;

macro_rules! testcase {
    ($name:ident) => {
        #[test]
        fn $name() {
            let name = stringify!($name);
            let testcase = load_testcase(
                &Path::new(file!())
                    .parent()
                    .unwrap()
                    .join("resources/testcases")
                    .join(name),
            )
            .unwrap();
            testcase.run().unwrap();
            testcase.check();
        }
    };
}

testcase!(test01);
testcase!(test02);
testcase!(test03);
testcase!(test04);
testcase!(test05);
testcase!(test06);
testcase!(test07);
testcase!(test08);
testcase!(test09);
testcase!(test10);
testcase!(test11);
testcase!(test12);
testcase!(test13);
testcase!(test14);
testcase!(test15);
testcase!(test16);
testcase!(test17);
testcase!(test18);
testcase!(test19);
testcase!(test20);
// skip the next test because this insertion cannot currently be resolved properly
// TODO find a way to fix this.
// testcase!(test21);
testcase!(test22);
testcase!(test23);
testcase!(test24);
testcase!(test25);
testcase!(test26);
testcase!(test27);
testcase!(test28);
testcase!(test29);
testcase!(test30);
testcase!(test31);
testcase!(test32);
testcase!(test33);
testcase!(test34);
testcase!(test35);
testcase!(test36);
testcase!(test37);
testcase!(test38);
testcase!(test39);
testcase!(test40);
testcase!(test41);
testcase!(test42);
testcase!(test43);
testcase!(test44);
testcase!(test45);
testcase!(test46);
testcase!(test47);
testcase!(test48);
testcase!(test49);
testcase!(test50);
testcase!(test51);
testcase!(pattern_too_long);
testcase!(test_wgbs01);
testcase!(test_long_pattern);
testcase!(test_contig_universe);

fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}

fn cleanup_file(f: &str) {
    if Path::new(f).exists() {
        fs::remove_file(f).unwrap();
    }
}

fn control_fdr(test: &str, event_str: &str, alpha: f64) {
    let basedir = basedir(test);
    let output = format!("{}/calls.filtered.bcf", basedir);
    cleanup_file(&output);
    varlociraptor::filtration::fdr::control_fdr(
        &format!("{}/calls.matched.bcf", basedir),
        Some(&output),
        &[varlociraptor::SimpleEvent {
            name: event_str.to_owned(),
        }],
        &varlociraptor::variants::model::VariantType::Deletion(Some(1..30)),
        LogProb::from(Prob(alpha)),
    )
    .unwrap();
}

fn assert_call_number(test: &str, expected_calls: usize) {
    let basedir = basedir(test);

    let mut reader = bcf::Reader::from_path(format!("{}/calls.filtered.bcf", basedir)).unwrap();

    let calls = reader.records().map(|r| r.unwrap()).collect_vec();
    // allow one more or less, in order to be robust to numeric fluctuations
    assert!(
        (calls.len() as i32 - expected_calls as i32).abs() <= 1,
        "unexpected number of calls ({} vs {})",
        calls.len(),
        expected_calls
    );
}

#[test]
fn test_fdr_control1() {
    control_fdr("test_fdr_ev_1", "SOMATIC", 0.05);
    //assert_call_number("test_fdr_ev_1", 974);
}

#[test]
fn test_fdr_control2() {
    control_fdr("test_fdr_ev_2", "SOMATIC", 0.05);
    assert_call_number("test_fdr_ev_2", 985);
}

/// same test, but low alpha
#[test]
fn test_fdr_control3() {
    control_fdr("test_fdr_ev_3", "ABSENT", 0.001);
    assert_call_number("test_fdr_ev_3", 0);
}

#[test]
fn test_fdr_control4() {
    control_fdr("test_fdr_ev_4", "SOMATIC_TUMOR", 0.05);
    assert_call_number("test_fdr_ev_4", 0);
}

// TODO enable this test again once https://github.com/samtools/bcftools/issues/874 is truly fixed upstream
// Then, also encode SVLEN as negative again for deletions.
//#[test]
fn test_fdr_control5() {
    control_fdr("test_fdr_control_out_of_bounds", "PRESENT", 0.05);
}

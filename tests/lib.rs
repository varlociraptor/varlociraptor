#[macro_use]
extern crate paste;
#[macro_use]
extern crate lazy_static;

use anyhow::Context;
use anyhow::Result;
use std::sync::Mutex;
use std::{fs, path::Path, path::PathBuf};

use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use rust_htslib::bcf::{self, Read, Reader};
use varlociraptor::{testcase, testcase_should_panic};

testcase!(test01, exact, fast);
testcase!(test02, exact, fast);
testcase!(test03, exact, fast);
testcase!(test04, exact, fast);
testcase!(test05, exact, fast);
testcase!(test06, exact, fast);
testcase!(test07, exact, fast);
testcase!(test08, exact, fast);
testcase!(test09, exact, fast);
testcase!(test10, exact, fast);
testcase!(test11, exact, fast);
testcase!(test12, exact, fast);
testcase!(test13, exact, fast);
testcase!(test14, exact, fast);
testcase!(test15, exact, fast);
testcase!(test16, exact, fast);
testcase!(test17, exact, fast);
testcase!(test18, exact, fast);
testcase!(test19, exact, fast);
testcase!(test20, exact, fast);
// skip the next test because this insertion cannot currently be resolved properly
// TODO find a way to fix this.
// testcase!(test21, exact, fast);
testcase!(test22, exact, fast);
testcase!(test23, exact, fast);
testcase!(test24, exact, fast);
testcase!(test25, exact, fast);
testcase!(test26, exact, fast);
testcase!(test27, exact, fast);
testcase!(test28, exact, fast);
testcase!(test29, exact, fast);
testcase!(test30, exact, fast);
testcase!(test31, exact, fast);
testcase!(test32, exact, fast);
testcase!(test33, exact, fast);
testcase!(test34, exact, fast);
testcase!(test36, exact, fast);
testcase!(test37, exact, fast);
// Skip this test. It is most likely a strand bias artifact, which is correctly recognized.
// However, there are also very few reads with nonstandard orientation, which are on the other
// strand.
//testcase!(test38, exact, fast);
testcase!(test39, exact, fast);
testcase!(test40, exact, fast);
testcase!(test41, exact, fast);
testcase!(test42, exact, fast);
testcase!(test43, exact, fast);
// Fast mode fails here, because there is a read with two insertions against
// the alt allele. This is very unlikely to happen, but it happens here.
// In the exact mode, there are various paths around this alignment which rescue
// the alt allele probability. With fast mode, these are missed, making the
// probability artificially small. That leads to Varlociraptor evaluating the
// locus to be heterozygous although it is homozygous in reality.
testcase!(test44, exact);
testcase!(test45, exact, fast);

testcase!(test47, exact, fast);
testcase!(test48, exact, fast);
testcase!(test49, exact);
testcase!(test50, exact);
testcase!(test51, exact);
testcase!(test52, exact, fast);
testcase!(test53, exact, fast);
testcase!(test54, exact, fast);
testcase!(test55, exact, fast);
testcase!(test57, exact);
testcase!(test58, exact);
testcase!(test59, exact);
testcase!(test60, exact);
testcase!(test61, exact);
testcase!(test62, exact);
testcase!(test63, exact);
testcase!(test64, exact);
testcase!(test65, exact);
testcase!(test69, exact);
testcase!(test70, exact);
testcase!(test71, exact);
testcase!(test72, exact);
testcase!(test74, exact);
testcase!(test75, exact);
testcase!(test77, exact);
testcase!(test78, exact);
testcase!(test79, exact);

// Skip test_giab_01 for now as the call seems to be correct.
// TODO try to find out what is wrong in the GIAB callset at that location.
testcase!(test_giab_01, exact);
testcase!(test_giab_02, exact);
testcase!(test_giab_03, exact);
testcase!(test_giab_04, exact);

testcase!(test_giab_05, exact);

testcase!(test_giab_06, exact);
testcase!(test_giab_07, exact);
testcase!(test_giab_08, exact);
testcase!(test_giab_09, exact);
testcase!(test_giab_10, exact);
testcase!(test_giab_11, exact);
testcase!(test_giab_12, exact);
testcase!(test_giab_13, exact);
// Skip test_giab_14. It is just bad luck that the reads here look like a homopolymer artifact although the variant is real.
// See testcase.yaml for details.
//testcase!(test_giab_14, exact);
testcase!(test_giab_15, exact);
testcase!(test_giab_16, exact);
testcase!(test_giab_17, exact);
testcase!(test_giab_18, exact);
testcase!(test_giab_19, exact);
testcase!(test_giab_20, exact);
testcase!(test_giab_21, exact);
testcase!(test_giab_22, exact);
testcase!(test_giab_23, exact);
// Skip test_giab_24. It is simply an unlucky combination of homology artifacts that looks strong.
// At least our probability is weaker than Freebayes's already.
// testcase!(test_giab_24, exact);
testcase!(test_giab_25, exact);
testcase!(test_giab_26, exact);
testcase!(test_giab_27, exact);
// Skip test_giab_28. It is simply an unlucky combination of homology artifacts that looks strongly like an artifact.
//testcase!(test_giab_28, exact);
testcase!(test_giab_29, exact);
testcase!(test_giab_30, exact);
testcase!(test_giab_31, exact);
testcase!(test_giab_32, exact);
testcase!(test_giab_33, exact);
testcase!(test_giab_34, exact);
testcase!(test_giab_35, exact);
testcase!(test_mapq_meth, exact);

testcase!(test_pcr_homopolymer_error1, exact);
testcase!(test_pcr_homopolymer_error2, exact);
testcase!(test_pcr_homopolymer_error3, exact);

testcase!(test_mendelian_prior, exact);
testcase!(pattern_too_long, exact, fast);
testcase!(test_long_pattern, exact, fast);
testcase!(test_contig_universe, exact, fast);
testcase!(test_expressions, exact);
testcase!(omit_sb, exact);
testcase!(test_panel_overlap, exact);
testcase!(test_panel_unknown_orientation_bias, exact);
testcase!(issue_154, exact, fast);
testcase!(test_low_cov_vaf, exact);
testcase_should_panic!(test_overlapping_events, exact);

testcase!(test_l2fc, exact, fast);
testcase!(test_cmp, exact, fast);

testcase!(test_nanopore_01, homopolymer);
testcase!(test_nanopore_02, homopolymer);
testcase!(test_nanopore_03, homopolymer);
testcase!(test_nanopore_04, homopolymer);
testcase!(test_nanopore_05, homopolymer);

testcase!(test_haplotype_absent, exact);
testcase!(test_haplotype_present, exact);
testcase!(test_haplotype_singleton, exact);

testcase!(test_alt_locus_bias_01, exact);
testcase!(test_alt_locus_bias_02, exact);
testcase!(test_uzuner_fn_mnv, exact);
testcase!(test_uzuner_fp_mnv1, exact);

testcase!(test_prinz_af_scan, exact);
testcase!(test_prinz_call_meth_1, exact);
testcase!(test_prinz_call_meth_2, exact);
testcase!(test_prinz_pacbio_zero, exact);
testcase!(test_prinz_ceta_one_sample_emseq, exact);
testcase!(test_prinz_ceta_one_sample_untreated, exact);
testcase!(test_prinz_ceta_two_samples, exact);

testcase!(test_imprecise_fusion, exact);
testcase!(test_imprecise_fusion_absent, exact);

testcase!(test_uzuner_clonal_1, exact);
testcase!(test_uzuner_clonal_2, exact);
testcase!(test_uzuner_clonal_3, exact);
testcase!(test_uzuner_fp_snv_on_ins, exact);
testcase!(test_false_negative_indel_call, exact);
testcase!(test_hiv_vaf_higher_than_expected, exact);
testcase!(test_uzuner_only_N, exact);
testcase!(test_moelder_floatisnan, exact);
testcase!(test_alt_locus_mapq_only, exact);

fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}

fn cleanup_file(f: &str) {
    if Path::new(f).exists() {
        fs::remove_file(f).unwrap();
    }
}

fn control_fdr(
    test: &str,
    events: &[&str],
    alpha: f64,
    local: bool,
    smart: bool,
    smart_retain_artifacts: bool,
    vartype: Option<&varlociraptor::variants::model::VariantType>,
) {
    let basedir = basedir(test);
    let output = format!("{}/calls.filtered.bcf", basedir);
    cleanup_file(&output);
    let event_strs: Vec<varlociraptor::SimpleEvent> = events
        .iter()
        .map(|&event_str| varlociraptor::SimpleEvent {
            name: event_str.to_owned(),
        })
        .collect();
    varlociraptor::filtration::fdr::control_fdr(
        &format!("{}/calls.matched.bcf", basedir),
        Some(&output),
        &event_strs,
        vartype,
        LogProb::from(Prob(alpha)),
        local,
        smart,
        smart_retain_artifacts,
    )
    .unwrap();
}

fn assert_call_number(test: &str, expected_calls: usize) {
    let basedir = basedir(test);

    let mut reader = bcf::Reader::from_path(format!("{}/calls.filtered.bcf", basedir)).unwrap();

    let calls = reader.records().map(|r| r.unwrap()).collect_vec();

    let ok = if expected_calls > 50 {
        // allow one more or less, in order to be robust to numeric fluctuations
        (calls.len() as i32 - expected_calls as i32).abs() <= 1
    } else {
        calls.len() == expected_calls
    };

    assert!(
        ok,
        "unexpected number of calls ({} vs {})",
        calls.len(),
        expected_calls
    );
}

#[test]
fn test_fdr_control1() {
    control_fdr(
        "test_fdr_ev_1",
        &["SOMATIC"],
        0.05,
        false,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    //assert_call_number("test_fdr_ev_1", 974);
}

#[test]
fn test_fdr_control2() {
    control_fdr(
        "test_fdr_ev_2",
        &["SOMATIC"],
        0.05,
        false,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_ev_2", 985);
}

/// same test, but low alpha
#[test]
fn test_fdr_control3() {
    control_fdr(
        "test_fdr_ev_3",
        &["ABSENT"],
        0.001,
        false,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_ev_3", 0);
}

#[test]
fn test_fdr_control4() {
    control_fdr(
        "test_fdr_ev_4",
        &["SOMATIC_TUMOR"],
        0.05,
        false,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_ev_4", 0);
}

#[test]
fn test_fdr_control_local1() {
    control_fdr(
        "test_fdr_local1",
        &["SOMATIC"],
        0.05,
        true,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_local1", 0);
}

#[test]
fn test_fdr_control_local2() {
    control_fdr(
        "test_fdr_local2",
        &["SOMATIC"],
        0.25,
        true,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_local2", 1);
}

#[test]
fn test_fdr_control_local2_smart() {
    control_fdr(
        "test_fdr_local2_smart",
        &["SOMATIC"],
        0.08,
        true,
        true,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_local2_smart", 1);
}

#[test]
fn test_fdr_control_local2_smart_retain_artifacts() {
    control_fdr(
        "test_fdr_local2_smart",
        &["SOMATIC"],
        0.08,
        true,
        true,
        true,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_local2_smart", 1);
}

#[test]
fn test_fdr_control_local3() {
    control_fdr(
        "test_fdr_local3",
        &["GERMLINE", "SOMATIC_TUMOR_LOW"],
        0.05,
        true,
        false,
        false,
        None,
    );
    assert_call_number("test_fdr_local3", 0);
}

// TODO enable this test again once https://github.com/samtools/bcftools/issues/874 is truly fixed upstream
// Then, also encode SVLEN as negative again for deletions.
//#[test]
// fn test_fdr_control5() {
//     control_fdr(
//         "test_fdr_control_out_of_bounds",
//         &["PRESENT"],
//         0.05,
//         false,
//         Some(&varlociraptor::variants::model::VariantType::Deletion(
//             Some(1..30),
//         )),
//     );
// }

//####################################################################################################################################################
// Tests for methylation candidates
//####################################################################################################################################################

fn control_meth_candidates(test: &str) -> Result<()> {
    let basedir = basedir(test);
    let output = format!("{}/candidates.bcf", basedir);
    cleanup_file(&output);
    varlociraptor::candidates::methylation::find_candidates(
        PathBuf::from(format!("{}/genome.fasta", basedir)),
        vec![varlociraptor::candidates::methylation::MethylationMotif::CG],
        Some(PathBuf::from(output)),
    )
    .with_context(|| "error computing methylation candidates".to_string())?;
    Ok(())
}

fn assert_candidates_number(test: &str, expected_calls: usize) -> Result<()> {
    let basedir = basedir(test);

    let mut reader = Reader::from_path(format!("{}/candidates.bcf", basedir))
        .with_context(|| "error reading BCF file".to_string())?;
    let calls = reader.records().map(|record| record.unwrap()).collect_vec();

    let ok = calls.len() == expected_calls;

    assert!(
        ok,
        "unexpected number of calls ({} vs {})",
        calls.len(),
        expected_calls
    );
    Ok(())
}

#[test]
fn test_meth_candidates1() -> Result<()> {
    control_meth_candidates("test_meth_ev_1")?;
    assert_candidates_number("test_meth_ev_1", 6)?;
    Ok(())
}

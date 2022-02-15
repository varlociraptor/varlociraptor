pub(crate) mod version0;
pub(crate) mod version1;
pub(crate) mod version2;
pub(crate) mod version3;
pub(crate) mod version4;

pub(crate) use version0::TestcaseVersion0;
pub(crate) use version1::TestcaseVersion1;
pub(crate) use version2::TestcaseVersion2;
pub(crate) use version3::TestcaseVersion3;
pub(crate) use version4::TestcaseVersion4;

use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::str;

use anyhow::Result;
use bio::io::fasta;
use eval::Expr;
use itertools::Itertools;
use rust_htslib::bcf::Read as BCFRead;
use rust_htslib::{bam, bcf};
use serde_json;
use tempfile::{self, NamedTempFile};
use yaml_rust::{Yaml, YamlLoader};

use varlociraptor::cli::{run, CallKind, PreprocessKind, VariantCallMode, Varlociraptor};
use varlociraptor::testcase::Mode;
use varlociraptor::utils;

pub(crate) fn load_testcase(path: impl AsRef<Path>) -> Result<Box<dyn Testcase>> {
    let mut reader = File::open(path.as_ref().join("testcase.yaml"))?;
    let mut content2 = String::new();
    reader.read_to_string(&mut content2)?;
    let yaml = YamlLoader::load_from_str(&content2)?;
    Ok(match &yaml[0]["version"] {
        Yaml::BadValue => Box::new(TestcaseVersion0 {
            inner: yaml,
            path: path.as_ref().to_owned(),
        }),
        Yaml::String(version) if version == "1" => Box::new(TestcaseVersion1 {
            inner: yaml,
            path: path.as_ref().to_owned(),
        }),
        Yaml::String(version) if version == "2" => Box::new(TestcaseVersion2 {
            inner: yaml,
            path: path.as_ref().to_owned(),
        }),
        Yaml::String(version) if version == "3" => Box::new(TestcaseVersion3 {
            inner: yaml,
            path: path.as_ref().to_owned(),
        }),
        Yaml::String(version) if version == "4" => Box::new(TestcaseVersion4 {
            inner: yaml,
            path: path.as_ref().to_owned(),
        }),
        _ => panic!("unsupported testcase version"),
    })
}

pub(crate) trait Testcase {
    fn inner(&self) -> &[Yaml];

    fn path(&self) -> &PathBuf;

    /// Index of record in the candidates.vcf to evaluate expressions for.
    fn test_record_index(&self) -> usize {
        self.yaml()
            .as_hash()
            .unwrap()
            .get(&Yaml::String("record-index".to_owned()))
            .map(|value| {
                value
                    .as_i64()
                    .expect("Invalid record index, expected integer") as usize
            })
            .unwrap_or(0)
    }

    fn preprocess_options(&self, sample_name: &str) -> String {
        self.yaml()["samples"][sample_name]["options"]
            .as_str()
            .unwrap()
            .to_owned()
    }

    fn mode(&self) -> Mode {
        match self.yaml()["mode"].as_str().unwrap() {
            "Generic" => Mode::Generic,
            "TumorNormal" => Mode::TumorNormal,
            _ => panic!("unsupported mode"),
        }
    }

    fn omit_strand_bias(&self) -> bool {
        if self.yaml()["omit_strand_bias"].is_badvalue() {
            false
        } else {
            self.yaml()["omit_strand_bias"].as_bool().unwrap()
        }
    }

    fn omit_read_orientation_bias(&self) -> bool {
        if self.yaml()["omit_read_orientation_bias"].is_badvalue() {
            false
        } else {
            self.yaml()["omit_read_orientation_bias"].as_bool().unwrap()
        }
    }

    fn omit_read_position_bias(&self) -> bool {
        if self.yaml()["omit_read_position_bias"].is_badvalue() {
            false
        } else {
            self.yaml()["omit_read_position_bias"].as_bool().unwrap()
        }
    }

    fn omit_softclip_bias(&self) -> bool {
        if self.yaml()["omit_softclip_bias"].is_badvalue() {
            false
        } else {
            self.yaml()["omit_softclip_bias"].as_bool().unwrap()
        }
    }

    fn omit_homopolymer_artifact_detection(&self) -> bool {
        if self.yaml()["omit_homopolymer_artifact_detection"].is_badvalue() {
            false
        } else {
            self.yaml()["omit_homopolymer_artifact_detection"]
                .as_bool()
                .unwrap()
        }
    }

    fn yaml(&self) -> &Yaml {
        &self.inner()[0]
    }

    fn output(&self) -> PathBuf {
        self.path().join("calls.bcf")
    }

    fn candidates(&self) -> PathBuf {
        self.path().join(self.yaml()["candidate"].as_str().unwrap())
    }

    fn samples(&self) -> Vec<String> {
        self.yaml()["samples"]
            .as_hash()
            .unwrap()
            .keys()
            .map(|name| name.as_str().unwrap().to_owned())
            .collect_vec()
    }

    fn sample_preprocessed_path(
        &self,
        sample_name: &str,
        temp_preprocess: &tempfile::TempDir,
    ) -> PathBuf {
        let mut path = temp_preprocess.as_ref().join(sample_name);
        path.set_extension("bcf");

        path
    }

    fn sample_observations_path(&self, sample_name: &str) -> PathBuf {
        let mut path = self.path().join(sample_name);
        path.set_extension("obs.tsv");

        path
    }

    fn sample(&self, sample_name: &str) -> &Yaml {
        &self.yaml()["samples"][sample_name]
    }

    fn sample_bam(&self, sample_name: &str) -> PathBuf {
        self.path()
            .join(self.sample(sample_name)["path"].as_str().unwrap())
    }

    fn sample_alignment_properties(&self, sample_name: &str) -> String {
        self.sample(sample_name)["properties"]
            .as_str()
            .unwrap()
            .to_owned()
    }

    fn scenario(&self) -> Option<PathBuf> {
        self.yaml()["scenario"]
            .as_str()
            .map(|p| self.path().join(p))
    }

    fn purity(&self) -> Option<f64> {
        self.yaml()["purity"].as_f64()
    }

    fn run(&self, pairhmm_mode_override: &str) -> Result<()> {
        let temp_ref = self.reference()?;

        let temp_preprocess = tempfile::tempdir()?;

        // Step 1: preprocess all samples
        for sample_name in &self.samples() {
            let mut options = serde_json::from_str(&self.preprocess_options(sample_name))?;
            match &mut options {
                Varlociraptor::Preprocess {
                    kind:
                        PreprocessKind::Variants {
                            ref mut reference,
                            ref mut candidates,
                            ref mut output,
                            ref mut bam,
                            ref mut alignment_properties,
                            ref mut output_raw_observations,
                            ..
                        },
                } => {
                    // prepare test bam
                    let test_bam = self.sample_bam(sample_name);
                    bam::index::build(&test_bam, None, bam::index::Type::Bai, 1).unwrap();

                    // prepare alignment properties
                    let props =
                        self.alignment_properties(&self.sample_alignment_properties(sample_name))?;

                    // replace options
                    *bam = test_bam;
                    *reference = PathBuf::from((*temp_ref).as_ref());
                    *candidates = self.candidates();
                    *output = Some(self.sample_preprocessed_path(sample_name, &temp_preprocess));
                    *alignment_properties = Some(props.path().to_owned());
                    *output_raw_observations = Some(self.sample_observations_path(sample_name));

                    run(options)?;
                }
                _ => panic!("bug: unsupported options"),
            }
        }

        // Step 2: run calling
        match self.mode() {
            Mode::Generic => {
                let options = Varlociraptor::Call {
                    kind: CallKind::Variants {
                        testcase_locus: None,
                        testcase_prefix: None,
                        testcase_anonymous: true,
                        omit_strand_bias: self.omit_strand_bias(),
                        omit_read_orientation_bias: self.omit_read_orientation_bias(),
                        omit_read_position_bias: self.omit_read_position_bias(),
                        omit_softclip_bias: self.omit_softclip_bias(),
                        omit_homopolymer_artifact_detection: self
                            .omit_homopolymer_artifact_detection(),
                        output: Some(self.output()),
                        mode: VariantCallMode::Generic {
                            scenario: self.scenario().unwrap(),
                            sample_observations: self
                                .samples()
                                .iter()
                                .map(|sample_name| {
                                    format!(
                                        "{}={}",
                                        sample_name,
                                        self.sample_preprocessed_path(
                                            sample_name,
                                            &temp_preprocess
                                        )
                                        .to_str()
                                        .unwrap()
                                    )
                                })
                                .collect_vec(),
                        },
                        log_mode: "default".to_owned(),
                    },
                };

                Ok(run(options)?)
            }
            Mode::TumorNormal => {
                let options = Varlociraptor::Call {
                    kind: CallKind::Variants {
                        testcase_locus: None,
                        testcase_prefix: None,
                        testcase_anonymous: true,
                        omit_strand_bias: self.omit_strand_bias(),
                        omit_read_orientation_bias: self.omit_read_orientation_bias(),
                        omit_read_position_bias: self.omit_read_position_bias(),
                        omit_softclip_bias: self.omit_softclip_bias(),
                        omit_homopolymer_artifact_detection: self
                            .omit_homopolymer_artifact_detection(),
                        output: Some(self.output()),
                        mode: VariantCallMode::TumorNormal {
                            tumor_observations: self
                                .sample_preprocessed_path("tumor", &temp_preprocess),
                            normal_observations: self
                                .sample_preprocessed_path("normal", &temp_preprocess),
                            purity: self.purity().unwrap(),
                        },
                        log_mode: "default".to_owned(),
                    },
                };

                Ok(run(options)?)
            }
        }
    }

    fn check(&self) {
        let mut reader = bcf::Reader::from_path(self.output()).unwrap();
        let mut calls = reader.records().map(|r| r.unwrap()).collect_vec();

        let calls: Box<dyn Iterator<Item = &mut bcf::Record>> =
            if utils::is_bnd(&mut calls[0]).expect("bug: failed to check for breakend") {
                Box::new(calls.iter_mut())
            } else {
                Box::new(calls.iter_mut().skip(self.test_record_index()).take(1))
            };

        for call in calls {
            let afs = call.format(b"AF").float().unwrap();
            if let Some(exprs) = self.yaml()["expected"]["allelefreqs"].as_vec() {
                for expr in exprs.iter() {
                    let mut expr = Expr::new(expr.as_str().unwrap());

                    for (sample, af) in reader.header().samples().into_iter().zip(afs.iter()) {
                        expr = expr.value(str::from_utf8(sample).unwrap(), af[0]);
                    }
                    assert!(
                        expr.exec()
                            .map(|v| v.as_bool().unwrap_or(false))
                            .unwrap_or(false),
                        "{:?} did not return true",
                        expr
                    );
                }
            }

            if let Some(exprs) = self.yaml()["expected"]["posteriors"].as_vec() {
                for expr in exprs.iter() {
                    let mut expr = Expr::new(expr.as_str().unwrap());

                    for rec in reader.header().header_records() {
                        match rec {
                            bcf::HeaderRecord::Info { values, .. } => {
                                let id = values.get("ID").unwrap().clone();
                                if id.starts_with("PROB_") {
                                    if let Ok(Some(values)) = call.info(id.as_bytes()).float() {
                                        expr = expr.value(id.clone(), values[0])
                                    }
                                }
                            }
                            _ => (), // ignore other tags
                        }
                    }
                    assert!(
                        expr.exec()
                            .map(|v| v.as_bool().unwrap_or(false))
                            .unwrap_or(false),
                        "{:?} did not return true",
                        expr
                    );
                }
            }
        }
    }

    fn reference(&self) -> Result<Box<dyn AsRef<Path>>> {
        let ref_name = self.yaml()["reference"]["name"].as_str().unwrap();
        let ref_seq = self.yaml()["reference"]["seq"].as_str().unwrap();

        let mut tmp_ref = tempfile::Builder::new().suffix(".fasta").tempfile()?;
        {
            let mut writer = fasta::Writer::new(&mut tmp_ref);
            writer.write(ref_name, None, ref_seq.as_bytes())?;
        }
        self.index_reference(&tmp_ref);

        Ok(Box::new(tmp_ref))
    }

    fn index_reference(&self, path: &dyn AsRef<Path>) {
        Command::new("samtools")
            .args(&["faidx", path.as_ref().to_str().unwrap()])
            .status()
            .expect("failed to create fasta index");
    }

    fn alignment_properties(&self, properties: &str) -> Result<NamedTempFile> {
        let mut tmp_props = tempfile::Builder::new().suffix(".json").tempfile()?;
        tmp_props.as_file_mut().write_all(properties.as_bytes())?;

        Ok(tmp_props)
    }
}

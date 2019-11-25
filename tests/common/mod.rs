pub mod version0;
pub mod version1;

pub use version0::TestcaseVersion0;
pub use version1::TestcaseVersion1;

use std::error::Error;

use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::str;

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

pub fn load_testcase(path: impl AsRef<Path>) -> Result<Box<dyn Testcase>, Box<dyn Error>> {
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
        _ => panic!("unsupported testcase version"),
    })
}

pub trait Testcase {
    fn inner(&self) -> &[Yaml];

    fn path(&self) -> &PathBuf;

    fn preprocess_options(&self, sample_name: &str) -> String;

    fn mode(&self) -> Mode;

    fn yaml(&self) -> &Yaml {
        &self.inner()[0]
    }

    fn output(&self) -> PathBuf {
        self.path().join("calls.bcf")
    }

    fn candidates(&self) -> PathBuf {
        self.path().join(self.yaml()["candidate"].as_str().unwrap())
    }

    fn reference_name(&self) -> &str {
        self.yaml()["reference"]["name"].as_str().unwrap()
    }

    fn reference_seq(&self) -> &str {
        self.yaml()["reference"]["seq"].as_str().unwrap()
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
        path.set_extension(".bcf");

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

    fn run(&self) -> Result<(), Box<dyn Error>> {
        let temp_ref = self.reference(self.reference_name(), self.reference_seq())?;

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
                            ..
                        },
                } => {
                    // prepare test bam
                    let test_bam = self.sample_bam(sample_name);
                    bam::index::build(&test_bam, None, bam::index::Type::BAI, 1).unwrap();

                    // prepare alignment properties
                    let props =
                        self.alignment_properties(&self.sample_alignment_properties(sample_name))?;

                    // replace options
                    *bam = test_bam;
                    *reference = temp_ref.path().to_owned();
                    *candidates = Some(self.candidates());
                    *output = Some(self.sample_preprocessed_path(sample_name, &temp_preprocess));
                    *alignment_properties = Some(props.path().to_owned());

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
                    },
                };

                Ok(run(options)?)
            }
            Mode::TumorNormal => {
                let options = Varlociraptor::Call {
                    kind: CallKind::Variants {
                        testcase_locus: None,
                        testcase_prefix: None,
                        output: Some(self.output()),
                        mode: VariantCallMode::TumorNormal {
                            tumor_observations: self
                                .sample_preprocessed_path("tumor", &temp_preprocess),
                            normal_observations: self
                                .sample_preprocessed_path("normal", &temp_preprocess),
                            purity: self.purity().unwrap(),
                        },
                    },
                };

                Ok(run(options)?)
            }
        }

        // let options = Varlociraptor::Call {
        //     kind: CallKind::Variants {

        //     }
        // }

        // match &mut options {
        //     Varlociraptor::Call { ref mut kind } => match kind {
        //         CallKind::Variants {
        //             ref mut mode,
        //             ref mut reference,
        //             ref mut candidates,
        //             ref mut output,
        //             ref mut testcase_locus,
        //             ref mut testcase_prefix,
        //             ..
        //         } => {
        //             *reference = temp_ref.path().to_owned();
        //             *candidates = Some(self.path.join(self.yaml()["candidate"].as_str().unwrap()));
        //             *output = Some(self.output());
        //             *testcase_prefix = None;
        //             *testcase_locus = None;

        //             match mode {
        //                 VariantCallMode::Generic {
        //                     ref mut scenario,
        //                     ref mut bams,
        //                     ref mut alignment_properties,
        //                 } => {
        //                     *scenario = self.path.join(self.yaml()["scenario"].as_str().unwrap());
        //                     bams.clear();
        //                     alignment_properties.clear();
        //                     let mut temp_props = Vec::new();
        //                     for (sample_name, sample) in
        //                         self.yaml()["samples"].as_hash().unwrap().iter()
        //                     {
        //                         let sample_name = sample_name.as_str().unwrap();
        //                         let bam = self.path.join(sample["path"].as_str().unwrap());
        //                         bam::index::build(&bam, None, bam::index::Type::BAI, 1).unwrap();
        //                         bams.push(format!("{}={}", sample_name, bam.to_str().unwrap()));
        //                         let props = Self::alignment_properties(
        //                             sample["properties"].as_str().unwrap(),
        //                         )?;
        //                         alignment_properties.push(format!(
        //                             "{}={}",
        //                             sample_name,
        //                             props.path().to_str().unwrap()
        //                         ));
        //                         temp_props.push(props);
        //                     }
        //                     run(options)
        //                 }
        //                 VariantCallMode::TumorNormal {
        //                     ref mut tumor,
        //                     ref mut normal,
        //                     ref mut tumor_alignment_properties,
        //                     ref mut normal_alignment_properties,
        //                     ..
        //                 } => {
        //                     *tumor = self
        //                         .path
        //                         .join(self.yaml()["samples"]["tumor"]["path"].as_str().unwrap());
        //                     *normal = self
        //                         .path
        //                         .join(self.yaml()["samples"]["normal"]["path"].as_str().unwrap());

        //                     let temp_tumor_props = Self::alignment_properties(
        //                         self.yaml()["samples"]["tumor"]["properties"]
        //                             .as_str()
        //                             .unwrap(),
        //                     )?;
        //                     let temp_normal_props = Self::alignment_properties(
        //                         self.yaml()["samples"]["normal"]["properties"]
        //                             .as_str()
        //                             .unwrap(),
        //                     )?;
        //                     *tumor_alignment_properties = Some(temp_tumor_props.path().to_owned());
        //                     *normal_alignment_properties =
        //                         Some(temp_normal_props.path().to_owned());

        //                     bam::index::build(tumor, None, bam::index::Type::BAI, 1).unwrap();
        //                     bam::index::build(normal, None, bam::index::Type::BAI, 1).unwrap();

        //                     run(options)
        //                 }
        //             }
        //         }
        //         _ => panic!("unsupported subcommand"),
        //     },
        //     _ => panic!("unsupported subcommand"),
        // }
    }

    fn check(&self) {
        let mut reader = bcf::Reader::from_path(self.output()).unwrap();
        let mut calls = reader.records().map(|r| r.unwrap()).collect_vec();
        assert_eq!(calls.len(), 1, "unexpected number of calls");
        let mut call = calls.pop().unwrap();

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
                                let values = call.info(id.as_bytes()).float().unwrap().unwrap();
                                expr = expr.value(id.clone(), values[0])
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

    fn reference(&self, ref_name: &str, ref_seq: &str) -> Result<NamedTempFile, Box<dyn Error>> {
        let mut tmp_ref = tempfile::Builder::new().suffix(".fasta").tempfile()?;
        {
            let mut writer = fasta::Writer::new(&mut tmp_ref);
            writer.write(ref_name, None, ref_seq.as_bytes())?;
        }
        Command::new("samtools")
            .args(&["faidx", tmp_ref.path().to_str().unwrap()])
            .status()
            .expect("failed to create fasta index");

        Ok(tmp_ref)
    }

    fn alignment_properties(&self, properties: &str) -> Result<NamedTempFile, Box<dyn Error>> {
        let mut tmp_props = tempfile::Builder::new().suffix(".json").tempfile()?;
        tmp_props.as_file_mut().write_all(properties.as_bytes())?;

        Ok(tmp_props)
    }
}

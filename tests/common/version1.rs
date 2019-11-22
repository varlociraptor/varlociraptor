use std::path::PathBuf;
use std::str;

use yaml_rust::Yaml;

use varlociraptor::testcase::Mode;
use crate::common::Testcase;

#[derive(Debug)]
pub struct TestcaseVersion1 {
    pub inner: Vec<Yaml>,
    pub path: PathBuf,
}

impl Testcase for TestcaseVersion1 {
    fn inner(&self) -> &[Yaml] {
        &self.inner
    }

    fn path(&self) -> &PathBuf {
        &self.path
    }

    fn mode(&self) -> Mode {
        match self.yaml()["mode"].as_str().unwrap() {
            "Generic" => Mode::Generic,
            "TumorNormal" => Mode::TumorNormal,
            _ => panic!("unsupported mode")
        }
    }

    fn preprocess_options(&self, sample_name: &str) -> String {
        self.yaml()["samples"][sample_name]["options"].as_str().unwrap().to_owned()
    }
}
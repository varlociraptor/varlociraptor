use std::path::PathBuf;

use yaml_rust::Yaml;

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
}

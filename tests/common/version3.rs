use std::path::{Path, PathBuf};

use anyhow::Result;
use serde_json::json;
use yaml_rust::Yaml;

use crate::common::Testcase;

#[derive(Debug)]
pub(crate) struct TestcaseVersion3 {
    pub(crate) inner: Vec<Yaml>,
    pub(crate) path: PathBuf,
}

impl Testcase for TestcaseVersion3 {
    fn inner(&self) -> &[Yaml] {
        &self.inner
    }

    fn path(&self) -> &PathBuf {
        &self.path
    }

    fn reference(&self) -> Result<Box<dyn AsRef<Path>>> {
        let reference_path = self.path.join(self.reference_path());
        self.index_reference(&reference_path);

        Ok(Box::new(reference_path.to_owned()))
    }

    fn preprocess_options(&self, sample_name: &str) -> String {
        let mut options: serde_json::Value = serde_json::from_str(
            self.yaml()["samples"][sample_name]["options"]
                .as_str()
                .unwrap(),
        )
        .unwrap();

        let variants = options["Preprocess"]["kind"].get_mut("Variants").unwrap();
        variants["candidates"] = json!("dummy.bcf");
        variants.as_object_mut().unwrap().remove("omit_snvs");
        variants.as_object_mut().unwrap().remove("omit_indels");

        options.to_string()
    }
}

impl TestcaseVersion3 {
    fn reference_path(&self) -> &Path {
        self.yaml()["reference"]["path"].as_str().unwrap().as_ref()
    }
}

use anyhow::Result;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::htslib::hts_idx_t;
use std::ffi;
use std::path::Path;

/// Get length, mapped and unmapped read counts for each contig.
pub(crate) fn idxstats(path: impl AsRef<Path>) -> Result<Vec<(String, (u64, u64, u64))>> {
    let mut reader = bam::Reader::from_path(path.as_ref())?;
    unsafe {
        let c_str = ffi::CString::new(path.as_ref().to_string_lossy().to_string())?;
        let htsfile = reader.htsfile();
        let idx = rust_htslib::htslib::sam_index_load(htsfile, c_str.as_ptr());
        Ok(_idxstats(&mut reader, idx))
    }
}

unsafe fn _idxstats<R: bam::Read>(
    bam: &mut R,
    idx: *mut hts_idx_t,
) -> Vec<(String, (u64, u64, u64))> {
    bam.header()
        .target_names()
        .iter()
        .enumerate()
        .map(|(tid, tname)| {
            let (mut mapped, mut unmapped) = (0, 0);
            rust_htslib::htslib::hts_idx_get_stat(idx, tid as i32, &mut mapped, &mut unmapped);
            let tname = std::str::from_utf8(tname).unwrap();
            let tlen = bam.header().target_len(tid as u32).unwrap();
            (tname.to_string(), (tlen, mapped, unmapped))
        })
        .collect()
}

//! Functions to check some properties of files and directories.
use std::{fs, io, path::Path};

/// Checks if the given folder is empty.
pub fn is_folder_empty(path: impl AsRef<Path>) -> io::Result<bool> {
    Ok(fs::read_dir(path)?.take(1).count() == 0)
}

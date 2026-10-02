//! Functions to remove files and directories.
use std::{fs, io, path::Path};

/// Cleans up the contents (files and folders) of the given folder while keeping the folder itself.
///
/// It is useful if you don't want to loose the permissions set on the folder, or if you only have
/// enough permissions to manipulate with the contents of the given folder, but not the folder
/// itself.
///
/// It is a common pattern:
///
/// * ["How to remove only the content of directories?"](https://unix.stackexchange.com/questions/45950/how-to-remove-only-the-content-of-directories)
/// * ["How to delete the contents of a folder in Python?"](https://stackoverflow.com/questions/185936/how-to-delete-the-contents-of-a-folder-in-python)
pub fn cleanup_folder(folder_path: impl AsRef<Path>) -> io::Result<()> {
    for entry in fs::read_dir(&folder_path)? {
        let path = entry?.path();
        if path.is_dir() {
            fs::remove_dir_all(&path)?;
        } else {
            fs::remove_file(&path)?;
        }
    }
    Ok(())
}

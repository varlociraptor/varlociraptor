//! Functions to copy files and directories from one place to another.
use quick_error::ResultExt;
use std::{
    fs, io,
    path::{Path, PathBuf},
};

struct SourceDirectory<'a>(&'a Path);

struct ObtainEntryIn<'a>(&'a Path);

struct CreateDirectory<'a>(&'a Path);

quick_error! {
    #[derive(Debug)]
    pub enum Error {
        CreateDirectory(p: PathBuf, err: io::Error) {
            description("A directory could not be created")
            display("Failed to create directory '{}'", p.display())
            context(p: CreateDirectory<'a>, err: io::Error) -> (p.0.to_path_buf(), err)
            cause(err)
        }
        ObtainEntry(p: PathBuf, err: io::Error) {
            description("A directory entry could not be obtained")
            display("Failed to read directory entry of '{}'", p.display())
            context(p: ObtainEntryIn<'a>, err: io::Error) -> (p.0.to_path_buf(), err)
            cause(err)
        }
        ReadDirectory(p: PathBuf, err: io::Error) {
            description("A directory could not be read to obtain its entries")
            display("Failed to read directory '{}'", p.display())
            context(p: SourceDirectory<'a>, err: io::Error) -> (p.0.to_path_buf(), err)
            cause(err)
        }
        DestinationDirectoryExists(p: PathBuf) {
            description("Cannot copy directories into an existing destination directory")
            display("Destination directory '{}' did already exist", p.display())
        }
        Copy(src: PathBuf, dest: PathBuf, err: io::Error) {
            description("A file could not be copied to its destination")
            display("Failed to copy '{}' to '{}'", src.display(), dest.display())
            context(c: (&'a PathBuf, &'a PathBuf), err: io::Error) -> (c.0.clone(), c.1.clone(), err)
            cause(err)
        }
    }
}

/// Return the computed destination directory, given a source directory.
pub fn destination_directory(
    source_dir: impl AsRef<Path>,
    destination_dir: impl AsRef<Path>,
) -> PathBuf {
    let source_dir = source_dir
        .as_ref()
        .canonicalize()
        .unwrap_or_else(|_| source_dir.as_ref().to_path_buf());
    destination_dir
        .as_ref()
        .join(source_dir.file_name().unwrap_or_else(|| "ROOT".as_ref()))
}

/// Copies the contents of the source directory to the given destination directory.
/// In `destination_dir`, a new subdirectory with the basename of the `source_dir` will be created.
/// It will not perform the copy operation if the effective destination directory does already exist.
///
/// The returned value will contain a copied directory's path.
pub fn copy_directory(
    source_dir: impl AsRef<Path>,
    destination_dir: impl AsRef<Path>,
) -> Result<PathBuf, Error> {
    let dest = destination_directory(source_dir.as_ref(), destination_dir);
    if dest.is_dir() {
        return Err(Error::DestinationDirectoryExists(dest));
    }

    // one possible implementation of walking a directory only visiting files
    fn visit_dirs(dir: &Path, dest: &Path) -> Result<(), Error> {
        for entry in fs::read_dir(dir).context(SourceDirectory(dir))? {
            let path = entry.context(ObtainEntryIn(dir))?.path();
            if path.is_dir() {
                visit_dirs(
                    &path,
                    &dest.join(path.file_name().expect("should always have filename here")),
                )?;
            } else {
                fs::create_dir_all(&dest).context(CreateDirectory(&dest))?;
                let dest = dest.join(&path.file_name().expect("should have filename here"));
                fs::copy(&path, &dest).context((&path, &dest))?;
            }
        }
        Ok(())
    }
    visit_dirs(source_dir.as_ref(), &dest)?;
    Ok(dest)
}

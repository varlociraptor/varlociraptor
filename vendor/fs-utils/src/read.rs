//! Functions to read from files.
use std::{
    self, fs,
    io::{self, Read},
    path::Path,
};

/// Reads the first N bytes from a file.
///
/// It is equivalent to `head -c limit` *nix utility.
pub fn head(path: impl AsRef<Path>, limit: usize) -> io::Result<Vec<u8>> {
    let file_size = fs::metadata(&path)?.len();
    let file_size: usize = if file_size <= std::usize::MAX as u64 {
        file_size as usize
    } else {
        std::usize::MAX
    };
    let (read_buffer_size, read_limit) = if file_size <= limit {
        (file_size, file_size)
    } else {
        (limit, limit)
    };
    let mut read_buffer = vec![0; read_buffer_size];
    fs::File::open(&path)?.read_exact(&mut read_buffer[..read_limit])?;
    Ok(read_buffer)
}

/// Reads the first `N` bytes from a file and return them as a string.
///
/// It assumes that the file is encoded with UTF-8, so any invalid UTF-8 sequences will be
/// replaced with `U+FFFD REPLACEMENT CHARACTER`, which looks like this: �, learn more
/// [here](https://doc.rust-lang.org/std/string/struct.string.html#method.from_utf8_lossy)).
///
/// It is equivalent to `head -c limit` *nix utility.
pub fn head_to_string(path: impl AsRef<Path>, limit: usize) -> io::Result<String> {
    Ok(String::from_utf8_lossy(&head(path, limit)?).into_owned())
}

/// Reads the first `N` bytes from a file and return them as a string. If the file size is greater
/// than `N` bytes, the truncation message will be put at the end of the String.
///
/// It assumes that the file is encoded with UTF-8, so any invalid UTF-8 sequences will be
/// replaced with `U+FFFD REPLACEMENT CHARACTER`, which looks like this: �, learn more
/// [here](https://doc.rust-lang.org/std/string/struct.string.html#method.from_utf8_lossy)).
pub fn head_to_string_with_message(
    path: impl AsRef<Path>,
    limit: usize,
    truncation_message: &str,
) -> io::Result<String> {
    let mut read_buffer = head(path, limit + 1)?;
    if read_buffer.len() > limit {
        read_buffer[(limit - truncation_message.len())..limit]
            .copy_from_slice(truncation_message.as_bytes());
    }
    read_buffer.truncate(limit);
    Ok(String::from_utf8_lossy(&read_buffer).into_owned())
}

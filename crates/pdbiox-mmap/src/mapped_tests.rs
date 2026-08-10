use super::*;
use std::fs::OpenOptions;
use std::io::{Seek, Write};

fn temporary_file(name: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!("pdbiox-map-{name}-{}", std::process::id()))
}

#[test]
fn safe_mapping_is_an_owned_snapshot() {
    let path = temporary_file("snapshot");
    let mut output = match OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(&path)
    {
        Ok(file) => file,
        Err(error) => panic!("create failed: {error}"),
    };
    if let Err(error) = output.write_all(b"mapped bytes") {
        panic!("write failed: {error}")
    }
    let mapped = match MappedFile::new(&output) {
        Ok(mapped) => mapped,
        Err(error) => panic!("snapshot failed: {error}"),
    };

    if let Err(error) = output.set_len(0) {
        panic!("truncate failed: {error}")
    }
    if let Err(error) = output.rewind() {
        panic!("rewind failed: {error}")
    }
    if let Err(error) = output.write_all(b"changed") {
        panic!("rewrite failed: {error}")
    }

    assert_eq!(mapped.as_bytes(), b"mapped bytes");
    let _removed = std::fs::remove_file(path);
}

#[test]
fn empty_file_has_an_empty_snapshot() {
    let path = temporary_file("empty");
    let input = match File::create(&path) {
        Ok(file) => file,
        Err(error) => panic!("create failed: {error}"),
    };
    let mapped = match MappedFile::new(&input) {
        Ok(mapped) => mapped,
        Err(error) => panic!("snapshot failed: {error}"),
    };
    assert!(mapped.is_empty());
    let _removed = std::fs::remove_file(path);
}

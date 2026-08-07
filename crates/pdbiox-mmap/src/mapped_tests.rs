use super::*;
use std::io::Write;

#[test]
fn a_mapping_borrows_file_pages_as_read_only_bytes() {
    let path = std::env::temp_dir().join(format!("pdbiox-map-{}", std::process::id()));
    let mut output = match File::create(&path) {
        Ok(file) => file,
        Err(error) => panic!("create failed: {error}"),
    };
    if let Err(error) = output.write_all(b"mapped bytes") {
        panic!("write failed: {error}")
    }
    drop(output);

    let input = match File::open(&path) {
        Ok(file) => file,
        Err(error) => panic!("open failed: {error}"),
    };
    let mapped = match MappedFile::new(&input) {
        Ok(mapped) => mapped,
        Err(error) => panic!("map failed: {error}"),
    };
    assert_eq!(mapped.as_bytes(), b"mapped bytes");
    let _removed = std::fs::remove_file(path);
}

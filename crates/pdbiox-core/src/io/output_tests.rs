use super::*;

#[test]
fn plain_output_writes_exactly_the_supplied_bytes() {
    let path = std::env::temp_dir().join(format!("pdbiox-output-{}.cif", std::process::id()));
    let bytes = b"data_test\n";
    if let Err(error) = write_output(&path, bytes) {
        panic!("write failed: {error}")
    }
    let actual = std::fs::read(&path).unwrap_or_default();
    let _removed = std::fs::remove_file(path);
    assert_eq!(actual, bytes);
}

#[cfg(feature = "gzip")]
#[test]
fn gzip_output_is_byte_deterministic() {
    let path = Path::new("test.cif.gz");
    let first = encode(path, b"data_test\n").unwrap_or_default();
    let second = encode(path, b"data_test\n").unwrap_or_default();
    assert_eq!(first, second);
}

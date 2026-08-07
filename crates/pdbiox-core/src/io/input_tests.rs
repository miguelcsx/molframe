use super::*;

#[test]
fn plain_bytes_are_returned_as_they_were_given() {
    let input = InputBuffer::from_bytes(b"data_test\n".to_vec());
    assert_eq!(input.as_bytes(), b"data_test\n");
    assert_eq!(input.len(), 10);
    assert!(!input.is_empty());
    assert_eq!(input.origin(), None);
}

#[test]
fn a_compression_container_is_recognised_from_its_leading_bytes() {
    assert_eq!(Compression::sniff(&[0x1f, 0x8b, 0x08]), Compression::Gzip);
    assert_eq!(
        Compression::sniff(&[0x28, 0xb5, 0x2f, 0xfd]),
        Compression::Zstd
    );
    assert_eq!(Compression::sniff(b"data_x"), Compression::None);
    assert_eq!(Compression::sniff(&[]), Compression::None);
}

#[test]
fn an_input_can_name_where_it_came_from_for_a_diagnostic_to_quote() {
    let input = InputBuffer::from_bytes(Vec::new()).with_origin("entry.cif");
    assert_eq!(input.origin(), Some("entry.cif"));
    assert!(input.is_empty());
}

#[test]
fn a_buffer_that_expands_beyond_the_limit_is_refused_rather_than_obeyed() {
    let limits = Limits {
        decompressed_bytes: 100,
        ..Limits::default()
    };
    assert!(check_expansion(101, 10, limits).is_err());
    assert!(check_expansion(99, 10, limits).is_ok());
}

#[test]
fn a_buffer_that_expands_faster_than_the_ratio_allows_is_refused() {
    let limits = Limits {
        compression_ratio: 10,
        ..Limits::default()
    };
    assert!(
        check_expansion(100, 10, limits).is_ok(),
        "tenfold is at the limit"
    );
    assert!(check_expansion(2_000, 10, limits).is_err());
}

#[test]
fn opening_a_file_that_is_not_there_reports_a_resource_finding() {
    let error = InputBuffer::open("no/such/entry.cif", Limits::default());
    assert_eq!(error.err().map(|finding| finding.code()), Some(Code::E1901));
}

#[test]
fn a_stream_is_read_and_bounded_like_a_file() {
    let input = InputBuffer::from_reader(std::io::Cursor::new(b"data_stream\n"), Limits::default());
    assert_eq!(
        input.as_ref().map(InputBuffer::as_bytes),
        Ok(b"data_stream\n".as_slice())
    );

    let limits = Limits {
        decompressed_bytes: 4,
        ..Limits::default()
    };
    assert!(InputBuffer::from_reader(std::io::Cursor::new(b"12345"), limits).is_err());
}

#[cfg(feature = "mmap")]
#[test]
fn a_large_plain_local_file_is_parsed_directly_from_a_mapping() {
    use std::io::Write as _;

    let path = std::env::temp_dir().join(format!("pdbiox-input-{}", std::process::id()));
    let mut file = match File::create(&path) {
        Ok(file) => file,
        Err(error) => panic!("create failed: {error}"),
    };
    let mut bytes = vec![b' '; MMAP_MIN_BYTES as usize];
    bytes[..10].copy_from_slice(b"data_test\n");
    if let Err(error) = file.write_all(&bytes) {
        panic!("write failed: {error}")
    }
    drop(file);

    let input = match InputBuffer::open(&path, Limits::default()) {
        Ok(input) => input,
        Err(finding) => panic!("open failed: {finding}"),
    };
    assert_eq!(input.kind(), InputKind::Mapped);
    assert!(input.as_bytes().starts_with(b"data_test\n"));
    let _removed = std::fs::remove_file(path);
}

#[cfg(feature = "mmap")]
#[test]
fn a_small_local_file_avoids_mapping_overhead() {
    use std::io::Write as _;

    let path = std::env::temp_dir().join(format!("pdbiox-small-input-{}", std::process::id()));
    let mut file = match File::create(&path) {
        Ok(file) => file,
        Err(error) => panic!("create failed: {error}"),
    };
    if let Err(error) = file.write_all(b"data_test\n") {
        panic!("write failed: {error}")
    }
    drop(file);

    let input = match InputBuffer::open(&path, Limits::default()) {
        Ok(input) => input,
        Err(finding) => panic!("open failed: {finding}"),
    };
    assert_eq!(input.kind(), InputKind::Owned);
    let _removed = std::fs::remove_file(path);
}

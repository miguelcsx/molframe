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

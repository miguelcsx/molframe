use super::*;

#[test]
fn decoding_plain_owned_input_preserves_its_payload_allocation() {
    let bytes = b"data_example\n".to_vec();
    let pointer = bytes.as_ptr();
    let input = InputBuffer::from_encoded_bytes(bytes, Limits::default()).unwrap();
    assert_eq!(input.as_bytes().as_ptr(), pointer);
}

#[test]
fn owned_encoded_input_obeys_the_absolute_byte_ceiling() {
    let limits = Limits {
        decompressed_bytes: 3,
        ..Limits::default()
    };
    assert!(InputBuffer::from_encoded_bytes(vec![0; 4], limits).is_err());
}

#[cfg(feature = "gzip")]
fn gzip(bytes: &[u8]) -> Vec<u8> {
    use std::io::Write;
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

#[cfg(feature = "gzip")]
#[test]
fn encoded_gzip_matches_stream_decoding_and_obeys_expansion_limits() {
    let expected = b"data_example\n".repeat(1000);
    let bytes = gzip(&expected);
    let streamed = InputBuffer::from_reader(bytes.as_slice(), Limits::default()).unwrap();
    let input = InputBuffer::from_encoded_bytes(bytes.clone(), Limits::default()).unwrap();
    assert_eq!(input.as_bytes(), expected);
    assert_eq!(input.as_bytes(), streamed.as_bytes());
    let limited = Limits {
        compression_ratio: 1,
        ..Limits::default()
    };
    assert!(InputBuffer::from_encoded_bytes(bytes, limited).is_err());
    assert!(InputBuffer::from_encoded_bytes(vec![0x1f, 0x8b, 0x08], Limits::default()).is_err());
}

#[cfg(feature = "zstd")]
#[test]
fn encoded_zstandard_matches_stream_decoding() {
    let expected = b"data_example\n".repeat(100);
    let bytes = zstd::stream::encode_all(expected.as_slice(), 1).unwrap();
    let input = InputBuffer::from_encoded_bytes(bytes, Limits::default()).unwrap();
    assert_eq!(input.as_bytes(), expected);
}

#[cfg(not(feature = "gzip"))]
#[test]
fn gzip_without_a_decoder_returns_an_error() {
    assert!(InputBuffer::from_encoded_bytes(vec![0x1f, 0x8b, 0x08], Limits::default()).is_err());
}

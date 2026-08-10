use super::*;

#[test]
fn all_integer_codec_families_decode() {
    assert_eq!(
        decode_i32(&blob(4, 2, 0, &i32s(&[4, -2]))).expect("codec 4"),
        [4, -2]
    );
    assert_eq!(
        decode_i32(&blob(7, 3, 0, &i32s(&[8, 2, 9, 1]))).expect("codec 7"),
        [8, 8, 9]
    );
    assert_eq!(
        decode_i32(&blob(8, 3, 0, &i32s(&[2, 3]))).expect("codec 8"),
        [2, 4, 6]
    );
    assert_eq!(
        decode_i32(&blob(14, 2, 0, &i16s(&[i16::MAX, 2, -4]))).expect("codec 14"),
        [32769, -4]
    );
    assert_eq!(
        decode_i32(&blob(15, 2, 0, &[127, 2, 252])).expect("codec 15"),
        [129, -4]
    );
}

#[test]
fn all_float_codec_families_decode() {
    let raw = [
        1.5f32.to_bits().to_be_bytes(),
        (-2.0f32).to_bits().to_be_bytes(),
    ]
    .concat();
    assert_eq!(
        decode_f32(&blob(1, 2, 0, &raw)).expect("codec 1"),
        [1.5, -2.0]
    );
    assert_eq!(
        decode_f32(&blob(9, 2, 10, &i32s(&[15, 2]))).expect("codec 9"),
        [1.5, 1.5]
    );
    assert_eq!(
        decode_f32(&blob(10, 2, 10, &i16s(&[10, 5]))).expect("codec 10"),
        [1.0, 1.5]
    );
    assert_eq!(
        decode_f32(&blob(11, 2, 10, &i16s(&[10, 15]))).expect("codec 11"),
        [1.0, 1.5]
    );
    assert_eq!(
        decode_f32(&blob(12, 2, 10, &i16s(&[10, 15]))).expect("codec 12"),
        [1.0, 1.5]
    );
    assert_eq!(
        decode_f32(&blob(13, 2, 10, &[10, 15])).expect("codec 13"),
        [1.0, 1.5]
    );
}

#[test]
fn encoded_output_round_trips_and_bad_lengths_refuse() {
    assert_eq!(
        decode_i32(&encode_i32(&[1, -7]).expect("encoding")).expect("integers"),
        [1, -7]
    );
    assert_eq!(
        decode_f32(&encode_f32(&[1.25, -0.5]).expect("encoding")).expect("floats"),
        [1.25, -0.5]
    );
    assert_eq!(
        decode_chars(&encode_chars(&[0, b'A']).expect("encoding")).expect("chars"),
        [0, b'A']
    );
    let strings = vec!["A".to_string(), "BC".to_string()];
    assert_eq!(
        decode_strings(&encode_strings(&strings, 4).expect("strings")).expect("decode"),
        strings
    );
    assert!(decode_i32(&blob(4, 2, 0, &i32s(&[1]))).is_err());
}

fn blob(kind: i32, length: usize, parameter: i32, payload: &[u8]) -> Vec<u8> {
    let mut bytes = header(kind, length, parameter).expect("header");
    bytes.extend_from_slice(payload);
    bytes
}

fn i32s(values: &[i32]) -> Vec<u8> {
    values
        .iter()
        .flat_map(|value| value.to_be_bytes())
        .collect()
}

fn i16s(values: &[i16]) -> Vec<u8> {
    values
        .iter()
        .flat_map(|value| value.to_be_bytes())
        .collect()
}

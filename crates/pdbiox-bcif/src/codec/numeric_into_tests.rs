use super::*;
use crate::{Decoded, decode, encode_floats, encode_interval};

fn assert_matches_allocating_decode(encoded: &EncodedData, expected: &[f64]) {
    let mut direct = vec![0.0_f32; expected.len()];
    let written = decode_f32_into(encoded, &mut direct).expect("direct numeric decode");
    assert_eq!(written, expected.len());
    let Decoded::Floats(reference) = decode(encoded).expect("reference decode") else {
        panic!("expected floating-point values");
    };
    for ((actual, reference), expected) in direct.iter().zip(reference).zip(expected) {
        assert!((f64::from(*actual) - reference).abs() <= f64::from(f32::EPSILON));
        assert!((f64::from(*actual) - expected).abs() <= f64::from(f32::EPSILON));
    }
}

#[test]
fn fixed_point_delta_or_packing_decodes_without_a_staging_column() {
    let values: Vec<f64> = (0_i16..300).map(|index| f64::from(index) * 0.125).collect();
    let encoded = encode_floats(&values).expect("encodable values");
    assert_matches_allocating_decode(&encoded, &values);
}

#[test]
fn interval_and_literal_float_chains_write_into_the_caller_buffer() {
    let values = [-3.0, -1.0, 0.5, 2.0, 9.0];
    let interval = encode_interval(&values, 4_096).expect("interval encoding");
    let Decoded::Floats(expected) = decode(&interval).expect("reference decode") else {
        panic!("expected floats");
    };
    assert_matches_allocating_decode(&interval, &expected);

    let literal =
        encode_floats(&[std::f64::consts::PI, std::f64::consts::E]).expect("literal encoding");
    assert_matches_allocating_decode(&literal, &[std::f64::consts::PI, std::f64::consts::E]);
}

#[test]
fn caller_output_length_is_enforced_before_any_hidden_growth() {
    let encoded = encode_floats(&[1.0, 2.0, 3.0]).expect("encoded values");
    let error = decode_f32_into(&encoded, &mut [0.0; 2]).expect_err("short output is refused");
    assert_eq!(error.code(), Code::E1401);
}

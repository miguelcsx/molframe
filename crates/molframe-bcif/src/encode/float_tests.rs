use super::*;
use crate::codec::{Decoded, decode};
use proptest::prelude::*;

#[test]
fn exact_float_encoding_round_trips_without_quantization() {
    let values = [1.25, -2.5, 3.0];
    let encoded = encode_floats(&values).expect("values encode");
    assert_eq!(
        decode(&encoded).expect("values decode"),
        Decoded::Floats(values.to_vec())
    );
}

proptest! {
    #[test]
    fn every_finite_float_column_round_trips_bits(values in proptest::collection::vec(any::<i32>(), 0..256)) {
        let values: Vec<f64> = values.into_iter().map(|value| f64::from(value) / 8.0).collect();
        let encoded = encode_floats(&values).expect("values encode");
        prop_assert_eq!(decode(&encoded).expect("values decode"), Decoded::Floats(values));
    }
}

#[test]
fn interval_encoding_keeps_endpoints_and_requested_bins() {
    let encoded = encode_interval(&[0.0, 0.5, 1.0], 3).expect("values encode");
    assert_eq!(
        decode(&encoded).expect("values decode"),
        Decoded::Floats(vec![0.0, 0.5, 1.0])
    );
}

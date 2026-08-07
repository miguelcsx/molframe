use super::*;
use crate::codec::{Decoded, decode};
use proptest::prelude::*;

#[test]
fn optimiser_round_trips_runs_deltas_and_large_values() {
    for values in [
        vec![7; 100],
        (0..100).map(i64::from).collect(),
        vec![-100_000, 0, 100_000, i32::MAX.into()],
    ] {
        let encoded = encode_integers(&values).expect("values encode");
        assert_eq!(
            decode(&encoded).expect("values decode"),
            Decoded::Integers(values)
        );
    }
}

#[test]
fn a_column_that_needs_signed_and_unsigned_thirty_two_bit_ranges_is_refused() {
    assert!(encode_integers(&[-1, i64::from(u32::MAX)]).is_err());
}

proptest! {
    #[test]
    fn every_integer_column_round_trips(values in proptest::collection::vec(any::<i32>(), 0..512)) {
        let values: Vec<i64> = values.into_iter().map(i64::from).collect();
        let encoded = encode_integers(&values).expect("values encode");
        prop_assert_eq!(decode(&encoded).expect("values decode"), Decoded::Integers(values));
    }
}

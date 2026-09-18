use super::*;
use crate::codec::{Decoded, decode};
use crate::encode::strategy::encoded_size;
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

#[test]
fn run_storage_is_built_only_when_its_byte_lower_bound_can_win() {
    let unique: Vec<i64> = (0..10_000).map(i64::from).collect();
    let encoding = Encoding::RunLength {
        src_type: DataType::Uint16,
        src_size: unique.len(),
    };
    assert!(
        runs_if_competitive(&unique, &encoding, unique.len())
            .expect("run count fits")
            .is_none()
    );

    let repeated = vec![7; 10_000];
    let runs = runs_if_competitive(&repeated, &encoding, repeated.len())
        .expect("run count fits")
        .expect("one run remains competitive");
    assert_eq!(runs, vec![7, 10_000]);
    assert_eq!(runs.capacity(), 2);
}

#[test]
fn integer_packing_writes_directly_into_exact_byte_capacity() {
    let values = [0, 255, 256];
    let (word_count, _) = packed_word_counts(&values, true).expect("word count fits");
    let bytes = pack_bytes(&values, 1, true, word_count).expect("values pack");
    assert_eq!(bytes, vec![0, 255, 0, 255, 1]);
    assert_eq!(bytes.capacity(), word_count);
}

#[test]
fn delta_storage_is_skipped_when_minimum_payload_cannot_win() {
    let values = vec![0; 10_000];
    let best = plain(&values, Range::of(&values)).expect("zero column has a plain encoding");
    let delta = Encoding::Delta {
        origin: 0,
        src_type: DataType::Uint8,
    };
    assert!(!can_outperform(&delta, values.len(), encoded_size(&best)));
}

proptest! {
    #[test]
    fn every_integer_column_round_trips(values in proptest::collection::vec(any::<i32>(), 0..512)) {
        let values: Vec<i64> = values.into_iter().map(i64::from).collect();
        let encoded = encode_integers(&values).expect("values encode");
        prop_assert_eq!(decode(&encoded).expect("values decode"), Decoded::Integers(values));
    }
}

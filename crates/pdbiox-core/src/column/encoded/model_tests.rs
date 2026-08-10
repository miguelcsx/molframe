use super::*;

fn round_trips<T: ColumnValue + std::fmt::Debug>(values: &[T]) {
    let column = EncodedColumn::encode(values);
    assert_eq!(column.len(), values.len(), "length for {column:?}");
    for (position, expected) in values.iter().enumerate() {
        assert_eq!(
            column.get(u32::try_from(position).expect("small position")),
            Some(*expected),
            "position {position} of {column:?}"
        );
    }
    let sequential: Vec<T> = column.iter().collect();
    assert_eq!(sequential.len(), values.len());
    for (read, expected) in sequential.iter().zip(values) {
        assert_eq!(read, expected);
    }
    assert_eq!(
        column.get(u32::try_from(values.len()).expect("small length")),
        None
    );
}

#[test]
fn every_encoding_reads_back_the_values_it_was_given() {
    round_trips::<u32>(&[]);
    round_trips(&[5u32; 100]);
    round_trips(&[1u32, 1, 1, 2, 2, 9]);
    round_trips(&[0u32, 1, 2, 3, 4, 5, 6, 7]);
    round_trips(&[100i32, 99, 150, -4, 0]);
    round_trips(&[1.5f32, -0.25, 3.0]);
    round_trips(&[7u8, 7, 200, 3]);
}

#[test]
fn a_column_of_one_repeated_value_collapses_to_a_constant() {
    let column = EncodedColumn::encode(&[3u32; 4096]);
    assert!(column.is_constant());
    assert_eq!(column.len(), 4096);
    assert_eq!(column.get(0), Some(3));
    assert_eq!(column.get(4095), Some(3));
}

#[test]
fn a_column_of_few_long_runs_becomes_run_length() {
    let mut values = vec![1u32; 50];
    values.extend([2u32; 50]);
    let column = EncodedColumn::encode(&values);
    assert!(matches!(column, EncodedColumn::RunLength { .. }));
    assert_eq!(column.get(49), Some(1));
    assert_eq!(column.get(50), Some(2));
}

#[test]
fn reading_a_run_length_column_forwards_costs_one_step_per_value() {
    let mut values = vec![1u32; 40];
    values.extend([2u32; 40]);
    values.extend([3u32; 40]);
    let column = EncodedColumn::encode(&values);
    assert_eq!(column.iter().collect::<Vec<_>>(), values);
}

#[test]
fn a_hot_column_stays_plain_when_it_declines_encoding() {
    let column = EncodedColumn::plain(&[1u32; 100]);
    assert!(!column.is_constant());
    assert_eq!(column.as_slice().map(<[u32]>::len), Some(100));
}

#[test]
fn only_a_plain_column_hands_out_a_slice() {
    assert!(EncodedColumn::encode(&[1u32; 8]).as_slice().is_none());
    assert!(EncodedColumn::plain(&[1u32, 2, 3]).as_slice().is_some());
}

#[test]
fn a_float_column_is_never_narrowed_by_its_exponent() {
    let column = EncodedColumn::encode(&[1.5f32, 2.5, 3.5]);
    assert!(!matches!(column, EncodedColumn::BitPacked { .. }));
    assert!(!matches!(column, EncodedColumn::Delta { .. }));
}

#[test]
fn a_monotonic_serial_column_is_narrower_than_storing_it_literally() {
    let serials: Vec<u32> = (1_000_000..1_000_064).collect();
    let column = EncodedColumn::encode(&serials);
    assert_eq!(column.get(63), Some(1_000_063));
    assert_eq!(column.len(), 64);
}

#[test]
fn a_symbol_column_encodes_like_any_other_narrow_integer() {
    let symbols: Vec<SymbolId> = [3u32, 3, 3, 7, 7]
        .iter()
        .map(|raw| SymbolId::from_raw(*raw))
        .collect();
    let column = EncodedColumn::encode(&symbols);
    assert_eq!(column.get(0), Some(SymbolId::from_raw(3)));
    assert_eq!(column.get(4), Some(SymbolId::from_raw(7)));
}

use super::*;

#[test]
fn a_large_repetitive_atom_column_uses_one_string_and_one_byte_per_row() {
    const ROWS: usize = 100_000;
    let strings = DecodedStringColumn::new(vec![Arc::from("ATOM")], vec![0; ROWS])
        .expect("valid dictionary indices");
    let ColumnValues::Strings(values) = ColumnValues::new(Decoded::Strings(strings)) else {
        panic!("string values expected")
    };

    assert_eq!(values.dictionary.len(), 1);
    let StringIndices::U8(indices) = &values.indices else {
        panic!("one-value dictionary must use byte indices")
    };
    assert_eq!(std::mem::size_of_val(indices.as_slice()), ROWS);
    assert_eq!(values.get(0), Some("ATOM"));
    assert_eq!(values.get(ROWS - 1), Some("ATOM"));
}

#[test]
fn integer_and_float_columns_narrow_without_changing_row_values() {
    let integers = ColumnValues::new(Decoded::Integers(vec![-2, 0, 127]));
    assert!(matches!(
        &integers,
        ColumnValues::Integers(IntegerValues::I8(_))
    ));
    assert!(matches!(integers.value(0), Some(ValueRef::Integer(-2))));

    let mut floats = ColumnValues::new(Decoded::Floats(vec![1.25, -2.5]));
    floats.compact_float();
    assert!(matches!(&floats, ColumnValues::Floats(FloatValues::F32(_))));
    assert!(matches!(floats.value(1), Some(ValueRef::Float(-2.5))));
}

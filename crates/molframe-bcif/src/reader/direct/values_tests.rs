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

#[test]
fn integers_pick_the_narrowest_width_that_holds_the_range() {
    let width = |values: Vec<i64>| match IntegerValues::new(values) {
        IntegerValues::I8(_) => 8,
        IntegerValues::I16(_) => 16,
        IntegerValues::I32(_) => 32,
        IntegerValues::I64(_) => 64,
    };
    assert_eq!(width(vec![-128, 127]), 8);
    assert_eq!(width(vec![0, 128]), 16);
    assert_eq!(width(vec![-32_769, 0]), 32);
    assert_eq!(width(vec![0, i64::from(i32::MAX) + 1]), 64);
    assert_eq!(width(Vec::new()), 8);
}

#[test]
fn string_indices_follow_the_dictionary_size_and_reject_stray_indices() {
    let column = |entries: usize, indices: Vec<u32>| {
        let dictionary: Vec<Arc<str>> = (0..entries).map(|n| Arc::from(n.to_string())).collect();
        let strings = DecodedStringColumn::new(dictionary, indices).expect("indices in range");
        let ColumnValues::Strings(values) = ColumnValues::new(Decoded::Strings(strings)) else {
            panic!("string values expected")
        };
        values
    };
    assert!(matches!(
        column(300, vec![299]).indices,
        StringIndices::U16(_)
    ));
    assert!(matches!(
        column(70_000, vec![69_999]).indices,
        StringIndices::U32(_)
    ));
    let small = column(10, vec![9, 0]);
    assert!(matches!(small.indices, StringIndices::U8(_)));
    assert_eq!(small.get(0), Some("9"));
    assert_eq!(small.get(1), Some("0"));
}

use super::*;
use arrow::ffi_stream::ArrowArrayStreamReader;

fn table() -> ColumnTable {
    match ColumnTable::new(vec![
        ("first".to_owned(), Column::U32(vec![1, 2, 3])),
        ("distance".to_owned(), Column::F32(vec![0.5, 1.5, 2.5])),
        ("score".to_owned(), Column::F64(vec![9.0, 8.0, 7.0])),
    ]) {
        Ok(table) => table,
        Err(error) => panic!("columns should bind: {error}"),
    }
}

#[test]
fn typed_columns_round_trip_through_a_c_stream() {
    let table = table();
    assert_eq!(table.len(), 3);
    let stream = match table.arrow_stream() {
        Ok(stream) => stream,
        Err(error) => panic!("stream creation failed: {error}"),
    };
    let reader = match ArrowArrayStreamReader::try_new(stream.into_ffi()) {
        Ok(reader) => reader,
        Err(error) => panic!("stream import failed: {error}"),
    };
    let batches: Vec<_> = reader
        .collect::<std::result::Result<_, _>>()
        .unwrap_or_default();
    assert_eq!(batches.len(), 1);
    let Some(distance) = batches[0]
        .column_by_name("distance")
        .and_then(|array| array.as_any().downcast_ref::<Float32Array>())
    else {
        panic!("distance column absent")
    };
    assert_eq!(distance.values(), &[0.5, 1.5, 2.5]);
    let schema = table.schema();
    assert_eq!(
        schema
            .field(0)
            .metadata()
            .get("molframe:export_cost")
            .map(String::as_str),
        Some("copy")
    );
}

#[test]
fn ragged_and_repeated_columns_are_refused() {
    let ragged = ColumnTable::new(vec![
        ("a".to_owned(), Column::U32(vec![1, 2])),
        ("b".to_owned(), Column::U32(vec![1])),
    ]);
    assert_eq!(
        ragged.err(),
        Some(ColumnTableError::LengthMismatch {
            name: "b".to_owned(),
            expected: 2,
            found: 1
        })
    );
    let repeated = ColumnTable::new(vec![
        ("a".to_owned(), Column::U8(vec![1])),
        ("a".to_owned(), Column::U8(vec![2])),
    ]);
    assert_eq!(
        repeated.err(),
        Some(ColumnTableError::DuplicateName {
            name: "a".to_owned()
        })
    );
}

#[test]
fn an_empty_table_has_a_schema_and_no_batches() {
    let empty = match ColumnTable::new(vec![("a".to_owned(), Column::I64(Vec::new()))]) {
        Ok(table) => table,
        Err(error) => panic!("an empty column binds: {error}"),
    };
    assert!(empty.is_empty());
    assert_eq!(empty.schema().fields().len(), 1);
    assert!(matches!(empty.record_batches(), Ok(batches) if batches.is_empty()));
}

use super::*;

fn table() -> ReflectionTable {
    let values = |items: &[i64]| {
        items
            .iter()
            .copied()
            .map(ReflectionValue::Integer)
            .collect()
    };
    ReflectionTable {
        title: "native".into(),
        cell: Some(UnitCell {
            lengths: [10.0, 11.0, 12.0],
            angles: [90.0, 91.0, 92.0],
        }),
        space_group_number: Some(1),
        space_group_name: Some("P 1".into()),
        columns: vec![
            ReflectionColumn {
                label: "H".into(),
                column_type: ReflectionColumnType::MillerIndex,
                values: values(&[0, 1]),
                dataset_id: 0,
                mtz_type: Some('H'),
            },
            ReflectionColumn {
                label: "K".into(),
                column_type: ReflectionColumnType::MillerIndex,
                values: values(&[0, 0]),
                dataset_id: 0,
                mtz_type: Some('H'),
            },
            ReflectionColumn {
                label: "L".into(),
                column_type: ReflectionColumnType::MillerIndex,
                values: values(&[1, 0]),
                dataset_id: 0,
                mtz_type: Some('H'),
            },
            ReflectionColumn {
                label: "FP".into(),
                column_type: ReflectionColumnType::Amplitude,
                values: vec![ReflectionValue::Real(10.5), ReflectionValue::Missing],
                dataset_id: 1,
                mtz_type: Some('F'),
            },
        ],
        datasets: vec![ReflectionDataset {
            id: 1,
            project: "proj".into(),
            crystal: "xtal".into(),
            name: "native".into(),
            wavelength: Some(1.0),
            cell: None,
        }],
        history: vec!["created by pdbiox".into()],
        symmetry_operations: vec!["X,Y,Z".into()],
        sort_order: [1, 2, 3, 0, 0],
        resolution_range: None,
        missing_value: None,
        extra_header_records: Vec::new(),
    }
}

#[test]
fn merged_mtz_round_trip_preserves_columns_datasets_and_history() {
    let expected = table();
    let bytes = write_mtz(&expected).expect("table should encode");
    let actual = read_mtz(&bytes).expect("MTZ should decode");
    assert_eq!(
        actual.miller_indices().expect("indices should exist"),
        [[0, 0, 1], [1, 0, 0]]
    );
    assert_eq!(actual.columns[3].values[1], ReflectionValue::Missing);
    assert_eq!(actual.datasets[0].name.as_ref(), "native");
    assert_eq!(actual.history, expected.history);
    assert_eq!(actual.symmetry_operations, expected.symmetry_operations);
}

#[test]
fn configured_gemmi_fixture_is_read_differentially() {
    let Some(path) = std::env::var_os("PDBIOX_MTZ_FIXTURE") else {
        return;
    };
    let bytes = std::fs::read(path).expect("configured Gemmi fixture should be readable");
    let table = read_mtz(&bytes).expect("Gemmi fixture should parse");
    assert_eq!(table.row_count(), 441);
    assert_eq!(table.columns.len(), 8);
    assert_eq!(
        table.miller_indices().expect("indices should exist")[0],
        [-5, 0, 1]
    );
}

#[test]
fn text_columns_are_refused_instead_of_silently_encoded() {
    let mut table = table();
    table.columns[3].values[0] = ReflectionValue::Text("observed".into());
    assert!(matches!(
        write_mtz(&table),
        Err(ReflectionError::Unsupported(_))
    ));
}

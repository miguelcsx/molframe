use super::*;
use molframe_cif::parse;
use molframe_core::io::InputBuffer;

const CIF: &str = "data_sf\n_cell.length_a 10\n_cell.length_b 11\n_cell.length_c 12\n_cell.angle_alpha 90\n_cell.angle_beta 91\n_cell.angle_gamma 92\n_symmetry.Int_Tables_number 4\n_symmetry.space_group_name_H-M 'P 1 21 1'\nloop_\n_refln.index_h\n_refln.index_k\n_refln.index_l\n_refln.F_meas_au\n_refln.F_meas_sigma_au\n_refln.status\n0 0 1 10.5 0.2 o\n1 0 0 ? . f\n";

fn document(text: &str) -> Document {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    parse(&input).expect("CIF should parse").0
}

#[test]
fn all_refln_columns_missing_kinds_and_metadata_are_lowered() {
    let table = lower_structure_factor_cif(&document(CIF)).expect("table should lower");
    assert_eq!(table.row_count(), 2);
    assert_eq!(
        table.miller_indices().expect("indices should exist")[1],
        [1, 0, 0]
    );
    assert_eq!(table.space_group_number, Some(4));
    assert_eq!(
        table.columns[3].column_type,
        ReflectionColumnType::Amplitude
    );
    assert_eq!(table.columns[3].values[1], ReflectionValue::Missing);
    assert_eq!(table.columns[4].values[1], ReflectionValue::Inapplicable);
    assert_eq!(
        table.columns[5].values[0],
        ReflectionValue::Text("o".into())
    );
}

#[test]
fn canonical_write_round_trips_the_reflection_model() {
    let mut expected = lower_structure_factor_cif(&document(CIF)).expect("table should lower");
    expected.symmetry_operations = vec!["X,Y,Z".into(), "-X,Y+1/2,-Z".into()];
    expected.datasets = vec![ReflectionDataset {
        id: 0,
        project: "".into(),
        crystal: "sf".into(),
        name: "sf".into(),
        wavelength: Some(1.234),
        cell: expected.cell,
    }];
    let text = write_structure_factor_cif(&expected).expect("table should write");
    let actual = lower_structure_factor_cif(&document(&text)).expect("written table should lower");
    assert_eq!(actual, expected);
}

#[test]
fn ragged_refln_loop_is_reported_by_the_cif_parser() {
    let input =
        InputBuffer::from_bytes(b"data_x\nloop_\n_refln.index_h\n_refln.index_k\n1\n".to_vec());
    let (document, diagnostics) = parse(&input).expect("tokenization should complete");
    assert!(!diagnostics.is_empty());
    assert!(document.first_block().is_some());
}

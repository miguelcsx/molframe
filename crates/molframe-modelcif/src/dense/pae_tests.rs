use super::PaeMatrix;
use crate::{PaeKind, PaePolicy, PaeStorage};

fn quarter_steps(size: usize) -> Vec<f32> {
    (0..size * size)
        .map(|index| (index % 127) as f32 * 0.25)
        .collect()
}

#[test]
fn stored_values_decode_to_exactly_the_supplied_values() {
    let values = quarter_steps(9);
    let matrix = PaeMatrix::from_values(9, &values, PaePolicy::new()).expect("matrix");
    assert!(matrix.is_exact());
    for (index, value) in values.iter().enumerate() {
        let (row, column) = (index / 9, index % 9);
        assert_eq!(matrix.get(row, column), Some(*value));
    }
}

#[test]
fn whole_number_matrices_are_one_byte_per_cell() {
    let values: Vec<f32> = (0..16).map(|index| f32::from(index as u8)).collect();
    let matrix = PaeMatrix::from_values(4, &values, PaePolicy::new()).expect("matrix");
    assert_eq!(matrix.kind(), PaeKind::U8);
    assert_eq!(matrix.byte_len(), 16);
    assert_eq!(matrix.as_u8().map(<[u8]>::len), Some(16));
    assert!(matrix.as_f32().is_none());
}

#[test]
fn forcing_one_byte_cells_reports_the_quantisation_error() {
    let values = [0.0, 31.75, 0.123_4, 20.0];
    let policy = PaePolicy::new().with_storage(PaeStorage::Quantised(PaeKind::U8));
    let matrix = PaeMatrix::from_values(2, &values, policy).expect("matrix");
    assert_eq!(matrix.kind(), PaeKind::U8);
    assert!(!matrix.is_exact());
    assert!(matrix.max_quantisation_error() <= 0.0625 + 1e-6);
}

#[test]
fn half_precision_is_chosen_when_it_is_the_only_exact_narrow_form() {
    let values = [2.0_f32.powi(-12), 1.5, 2.0, 3.0];
    let matrix = PaeMatrix::from_values(2, &values, PaePolicy::new()).expect("matrix");
    assert!(matrix.is_exact());
    assert_eq!(matrix.kind(), PaeKind::F16);
    assert_eq!(matrix.byte_len(), 8);
}

#[test]
fn an_oversized_matrix_is_refused_by_the_budget() {
    let values = quarter_steps(8);
    let policy = PaePolicy::new().with_max_bytes(10);
    assert!(PaeMatrix::from_values(8, &values, policy).is_err());
}

#[test]
fn rows_and_blocks_are_decoded_in_place() {
    let values = quarter_steps(4);
    let matrix = PaeMatrix::from_values(4, &values, PaePolicy::new()).expect("matrix");
    let mut row = [0.0; 4];
    matrix.row_into(2, &mut row).expect("row");
    assert_eq!(row, [values[8], values[9], values[10], values[11]]);
    let block = matrix.block_mean(0..2, 0..2).expect("block");
    let expected = f64::from(values[0] + values[1] + values[4] + values[5]) / 4.0;
    assert!((block - expected).abs() < 1e-9);
    assert!(matrix.block_mean(0..5, 0..1).is_none());
}

#[test]
fn setting_a_cell_outside_the_matrix_is_an_error() {
    let mut matrix = PaeMatrix::from_values(2, &[0.0; 4], PaePolicy::new()).expect("matrix");
    assert!(matrix.set(2, 0, 1.0).is_err());
    assert!(matrix.set(0, 0, f32::NAN).is_err());
}

fn matrix_from_f64(size: usize, values: &[f64]) -> PaeMatrix {
    let mut probe = crate::PrecisionProbe::new();
    for value in values {
        probe.observe_f64(*value);
    }
    let policy = PaePolicy::new();
    let layout = policy.layout(&probe, size).expect("layout");
    let mut matrix = PaeMatrix::zeroed(size, layout, policy.max_bytes).expect("matrix");
    for (index, value) in values.iter().enumerate() {
        matrix
            .set_f64(index / size, index % size, *value)
            .expect("cell");
    }
    matrix
}

#[test]
fn decimal_text_values_are_reproduced_exactly_in_binary64() {
    let values = [12.34, 0.07, 31.75, 0.0];
    let matrix = matrix_from_f64(2, &values);
    assert_eq!(matrix.kind(), PaeKind::U16);
    assert!(matrix.is_exact());
    for (index, value) in values.iter().enumerate() {
        assert_eq!(matrix.get_f64(index / 2, index % 2), Some(*value));
    }
}

#[test]
fn digits_beyond_single_precision_are_reported_as_an_error() {
    let values = [1.234_567_890_123, 2.0, 3.0, 4.0];
    let matrix = matrix_from_f64(2, &values);
    assert_eq!(matrix.kind(), PaeKind::F32);
    assert!(!matrix.is_exact());
    assert!(matrix.max_quantisation_error() < 1e-6);
}

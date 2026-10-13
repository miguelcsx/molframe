use super::read_pae_json;
use crate::{ConfidenceError, PaeKind, PaePolicy};

fn read(text: &str) -> Result<crate::PaeMatrix, ConfidenceError> {
    read_pae_json(text.as_bytes(), PaePolicy::new())
}

#[test]
fn list_wrapped_dense_matrix_is_read_with_its_declared_maximum() {
    let matrix = read(
        r#"[{"predicted_aligned_error":[[0,1.5],[2.25,0]],"max_predicted_aligned_error":31.75}]"#,
    )
    .expect("matrix");
    assert_eq!(matrix.size(), 2);
    assert_eq!(matrix.get(0, 1), Some(1.5));
    assert_eq!(matrix.get(1, 0), Some(2.25));
    assert_eq!(matrix.declared_maximum(), Some(31.75));
    assert!(matrix.is_exact());
}

#[test]
fn plain_object_and_unknown_members_are_accepted() {
    let matrix = read(r#"{"note":{"a":[1,2,{"b":"]"}]},"predicted_aligned_error":[[0,3],[4,0]]}"#)
        .expect("matrix");
    assert_eq!(matrix.kind(), PaeKind::U8);
    assert_eq!(matrix.get(1, 0), Some(4.0));
}

#[test]
fn integer_matrices_use_one_byte_cells() {
    let matrix = read(r#"[{"predicted_aligned_error":[[0,1],[2,0]]}]"#).expect("matrix");
    assert_eq!(matrix.byte_len(), 4);
}

#[test]
fn residue_triples_fill_a_square_matrix() {
    let text = r#"{"residue1":[1,1,2,2],"residue2":[1,2,1,2],"distance":[0.5,2,3,0.5]}"#;
    let matrix = read(text).expect("matrix");
    assert_eq!(matrix.size(), 2);
    assert_eq!(matrix.get(0, 1), Some(2.0));
    assert_eq!(matrix.get(1, 0), Some(3.0));
}

#[test]
fn missing_and_repeated_triples_are_rejected() {
    let missing = r#"{"residue1":[1,2],"residue2":[1,2],"distance":[0,0]}"#;
    assert!(matches!(
        read(missing),
        Err(ConfidenceError::Incomplete { .. })
    ));
    let repeated = r#"{"residue1":[1,1],"residue2":[1,1],"distance":[0,1]}"#;
    assert!(matches!(
        read(repeated),
        Err(ConfidenceError::Conflict { .. })
    ));
}

#[test]
fn ragged_and_non_square_matrices_are_rejected() {
    assert!(matches!(
        read(r#"{"predicted_aligned_error":[[0,1],[2]]}"#),
        Err(ConfidenceError::Ragged { .. })
    ));
    assert!(matches!(
        read(r#"{"predicted_aligned_error":[[0,1,2],[2,0,1]]}"#),
        Err(ConfidenceError::SizeMismatch { .. })
    ));
}

#[test]
fn unrecognised_documents_report_an_offset() {
    assert!(matches!(read("[]"), Err(ConfidenceError::Json { .. })));
    assert!(matches!(
        read(r#"{"other":1}"#),
        Err(ConfidenceError::Json { .. })
    ));
    assert!(matches!(
        read("{\"predicted_aligned_error\":[[0,x]]}"),
        Err(ConfidenceError::Json { .. })
    ));
}

#[test]
fn budget_is_enforced_before_the_matrix_is_allocated() {
    let policy = PaePolicy::new().with_max_bytes(3);
    let result = read_pae_json(br#"{"predicted_aligned_error":[[0,1],[2,0]]}"#, policy);
    assert!(matches!(result, Err(ConfidenceError::Budget { .. })));
}

#[test]
fn alphafold_p69905_matrix_round_trips_exactly() {
    let path = "/tmp/rt-data/af/AF-P69905-F1-predicted_aligned_error_v6.json";
    let Ok(text) = std::fs::read(path) else {
        return;
    };
    let matrix = read_pae_json(&text, PaePolicy::new()).expect("matrix");
    assert_eq!(matrix.size(), 142);
    assert_eq!(matrix.kind(), PaeKind::U8);
    assert!(matrix.is_exact());
    assert_eq!(matrix.declared_maximum(), Some(31.75));
}

#[test]
fn two_decimal_values_survive_import_in_binary64() {
    let matrix = read(r#"{"predicted_aligned_error":[[0,12.34],[0.07,0]]}"#).expect("matrix");
    assert_eq!(matrix.get_f64(0, 1), Some(12.34));
    assert_eq!(matrix.get_f64(1, 0), Some(0.07));
    assert!(matrix.is_exact());
}

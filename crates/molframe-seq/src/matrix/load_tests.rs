use super::load_matrix;
use crate::matrix::ScoreMatrix;
use molframe_chem::{TableError, TableRequest};

fn load(name: &str) -> crate::SubstitutionMatrix {
    let request = TableRequest::<ScoreMatrix>::Named(name.into());
    load_matrix(&request, molframe_bench::table_store()).expect("matrix loads")
}

#[test]
fn every_distributed_matrix_loads_with_its_own_name() {
    let mut names: Vec<String> = (30..=90).step_by(5).map(|n| format!("blosum{n}")).collect();
    names.extend(["identity".to_owned(), "nuc44".to_owned()]);
    for name in names {
        let matrix = load(&name);
        assert_eq!(matrix.identity().name, name);
    }
}

#[test]
fn nuc44_scores_matches_and_mismatches() {
    let matrix = load("nuc44");
    assert_eq!(matrix.get(b'A', b'A'), 5);
    assert_eq!(matrix.get(b'A', b'T'), -4);
}

#[test]
fn identity_matrix_rewards_only_exact_matches() {
    let matrix = load("identity");
    assert_eq!(matrix.get(b'W', b'W'), 1);
    assert_eq!(matrix.get(b'W', b'A'), -10000);
}

#[test]
fn a_pam_matrix_carries_its_provenance() {
    let matrix = load("pam250");
    assert!(matrix.identity().source.contains("ncbi"));
    assert!(!matrix.identity().version.is_empty());
}

#[test]
fn a_missing_matrix_name_is_an_absent_table() {
    let request = TableRequest::<ScoreMatrix>::Named("no-such-matrix".into());
    let error = load_matrix(&request, molframe_bench::table_store()).expect_err("absent");
    assert!(matches!(
        error,
        TableError::Absent { .. } | TableError::Missing { .. }
    ));
}

#[test]
fn a_malformed_file_is_reported_as_malformed() {
    let dir = std::env::temp_dir().join("molframe-seq-bad-matrix");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("bad.matrix");
    std::fs::write(
        &path,
        "format: 1\nname: bad\nversion: 1\nsource: s\nlicence: l\nA R\nA 1\n",
    )
    .expect("write");
    let request = TableRequest::<ScoreMatrix>::Named(path.to_string_lossy().into_owned());
    let error = load_matrix(&request, molframe_bench::table_store()).expect_err("malformed");
    assert!(matches!(error, TableError::Malformed { .. }));
}

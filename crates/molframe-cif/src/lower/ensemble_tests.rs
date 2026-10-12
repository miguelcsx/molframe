use super::*;
use crate::parse;
use molframe_core::io::InputBuffer;

fn atom_site(text: &str) -> crate::Category {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (document, _) = match parse(&input) {
        Ok(parsed) => parsed,
        Err(findings) => panic!("parse failed: {findings:?}"),
    };
    let Some(category) = document
        .first_block()
        .and_then(|block| block.category("atom_site"))
        .cloned()
    else {
        panic!("atom_site missing")
    };
    category
}

const HEADER: &str = "data_x\nloop_\n_atom_site.id\n_atom_site.label_atom_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\n";

#[test]
fn coordinate_changes_alone_are_dense() {
    let category = atom_site(&format!("{HEADER}1 N 0 0 0 1\n1 N 1 1 1 2\n#\n"));
    assert!(ragged_model_numbers(&category, false).is_empty());
}

#[test]
fn an_identity_change_is_ragged_even_when_atom_counts_match() {
    let category = atom_site(&format!("{HEADER}1 N 0 0 0 1\n1 O 1 1 1 2\n#\n"));
    assert_eq!(ragged_model_numbers(&category, false), [1, 2]);
}

#[test]
fn a_different_atom_count_is_ragged() {
    let category = atom_site(&format!(
        "{HEADER}1 N 0 0 0 1\n2 CA 1 0 0 1\n1 N 1 1 1 2\n#\n"
    ));
    assert_eq!(ragged_model_numbers(&category, false), [1, 2]);
}

#[test]
fn a_serial_continued_across_models_is_dense() {
    let category = atom_site(&format!(
        "{HEADER}1 N 0 0 0 1\n2 CA 1 0 0 1\n3 N 1 1 1 2\n4 CA 2 0 0 2\n#\n"
    ));
    assert!(ragged_model_numbers(&category, false).is_empty());
}

#[test]
fn a_continued_serial_does_not_hide_an_identity_change() {
    let category = atom_site(&format!(
        "{HEADER}1 N 0 0 0 1\n2 CA 1 0 0 1\n3 N 1 1 1 2\n4 CB 2 0 0 2\n#\n"
    ));
    assert_eq!(ragged_model_numbers(&category, false), [1, 2]);
}

#[test]
fn per_model_occupancy_and_b_factor_alone_stay_dense_and_are_flagged() {
    let header = "data_x\nloop_\n_atom_site.id\n_atom_site.label_atom_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_PDB_model_num\n";
    let differing = atom_site(&format!("{header}1 N 0 0 0 1 10 1\n1 N 1 1 1 0.5 20 2\n#\n"));
    assert!(ragged_model_numbers(&differing, false).is_empty());
    assert!(per_model_values_differ(&differing, false));
    assert!(!per_model_values_differ(&differing, true));
    let same = atom_site(&format!("{header}1 N 0 0 0 1 10 1\n1 N 1 1 1 1 10 2\n#\n"));
    assert!(!per_model_values_differ(&same, false));
}

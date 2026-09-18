use super::{EmptyQsPolicy, QsOptions, qs_score};
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(body: &str) -> Structure {
    let input = InputBuffer::from_bytes(format!("{HEADER}{body}").into_bytes());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

const NATIVE: &str = "\
ATOM 1 C CA GLY A 1 0 0 0\n\
ATOM 2 C CA GLY A 2 2 0 0\n\
ATOM 3 C CA GLY B 1 1 1 3\n\
ATOM 4 C CA GLY B 2 3 1 3\n";

#[test]
fn an_identical_interface_scores_one() {
    let native = structure(NATIVE);
    let model = structure(NATIVE);
    let score = qs_score(&model, &native, "A", "B", QsOptions::standard(5.0))
        .unwrap_or_else(|error| panic!("valid QS: {error}"));
    assert!((score - 1.0).abs() < 1e-12, "score {score}");
}

#[test]
fn a_pulled_apart_interface_scores_zero() {
    let native = structure(NATIVE);
    let model = structure(
        "ATOM 1 C CA GLY A 1 0 0 0\n\
ATOM 2 C CA GLY A 2 2 0 0\n\
ATOM 3 C CA GLY B 1 51 1 3\n\
ATOM 4 C CA GLY B 2 53 1 3\n",
    );
    let score = qs_score(&model, &native, "A", "B", QsOptions::standard(5.0))
        .unwrap_or_else(|error| panic!("valid QS: {error}"));
    assert!(score.abs() < 1e-12, "score {score}");
}

#[test]
fn empty_interface_and_cutoff_policies_are_explicit() {
    let native = structure(NATIVE);
    assert!(matches!(
        qs_score(
            &native,
            &native,
            "missing-a",
            "missing-b",
            QsOptions {
                contact_distance: 5.0,
                empty_policy: EmptyQsPolicy::Error,
            }
        ),
        Err(crate::CompareError::NoComparablePairs)
    ));
    assert!(matches!(
        qs_score(&native, &native, "A", "B", QsOptions::standard(f32::NAN)),
        Err(crate::CompareError::InvalidDistanceCutoff)
    ));
}

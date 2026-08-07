use super::*;

#[test]
fn coordinate_changes_alone_are_dense() {
    let text = "MODEL        4\nATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N\nENDMDL\nMODEL        8\nATOM      1  N   GLY A   1      2.000   2.000   2.000  1.00  0.00           N\nENDMDL\n";
    assert!(ragged_models(text, false).is_none());
}

#[test]
fn atom_identity_and_count_changes_are_ragged() {
    let identity = "MODEL        1\nATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N\nENDMDL\nMODEL        2\nATOM      1  O   GLY A   1      2.000   2.000   2.000  1.00  0.00           O\nENDMDL\n";
    let count = format!(
        "{identity}MODEL        3\nATOM      1  O   GLY A   1      3.000   3.000   3.000  1.00  0.00           O\nATOM      2  C   GLY A   1      4.000   4.000   4.000  1.00  0.00           C\nENDMDL\n"
    );
    assert_eq!(
        ragged_models(identity, false).map(|models| models.len()),
        Some(2)
    );
    assert_eq!(
        ragged_models(&count, false).map(|models| models.len()),
        Some(3)
    );
}

use super::*;

fn key(chain: u32, label: Option<i32>, auth: Option<i32>, ins: u32) -> ResidueKey {
    ResidueKey {
        model: 1,
        chain,
        label_seq: OptionalI32::from(label),
        auth_seq: OptionalI32::from(auth),
        ins_code: ins,
    }
}

#[test]
fn a_different_sequence_position_starts_a_new_residue() {
    let first = key(1, Some(1), Some(1), ResidueKey::ABSENT);
    let second = key(1, Some(2), Some(2), ResidueKey::ABSENT);
    assert_eq!(boundary(&first, &second, false), Boundary::New);
}

#[test]
fn a_different_chain_starts_a_new_residue() {
    let first = key(1, Some(1), Some(1), ResidueKey::ABSENT);
    let second = key(2, Some(1), Some(1), ResidueKey::ABSENT);
    assert_eq!(boundary(&first, &second, false), Boundary::New);
}

#[test]
fn an_insertion_code_is_part_of_identity_rather_than_an_annotation_on_it() {
    let plain = key(1, Some(163), Some(163), ResidueKey::ABSENT);
    let inserted = key(1, Some(163), Some(163), 7);
    assert_eq!(boundary(&plain, &inserted, false), Boundary::New);
    assert_ne!(plain, inserted);
}

#[test]
fn identical_annotations_continue_the_same_residue_by_default() {
    let first = key(1, Some(1), Some(1), ResidueKey::ABSENT);
    let second = key(1, Some(1), Some(1), ResidueKey::ABSENT);
    assert_eq!(boundary(&first, &second, false), Boundary::Same);
}

#[test]
fn a_repeated_atom_name_under_identical_annotations_falls_back_to_file_order() {
    let first = key(1, Some(1), Some(1), ResidueKey::ABSENT);
    let second = key(1, Some(1), Some(1), ResidueKey::ABSENT);
    assert_eq!(boundary(&first, &second, true), Boundary::NewByFileOrder);
}

#[test]
fn the_component_code_is_not_part_of_the_key_so_a_modelled_mutation_stays_one_residue() {
    // Two alternate locations of one residue, one carrying a different
    // component code. The key cannot see the code, so both land together.
    let conformer_a = key(1, Some(57), Some(57), ResidueKey::ABSENT);
    let conformer_b = key(1, Some(57), Some(57), ResidueKey::ABSENT);
    assert_eq!(boundary(&conformer_a, &conformer_b, false), Boundary::Same);
}

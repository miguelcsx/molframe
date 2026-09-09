use super::*;
use pdbiox_core::io::InputBuffer;
use pdbiox_core::structure::AtomRef;
use std::fmt::Write as _;

#[test]
fn common_records_do_not_retain_atom_lines() {
    let text = "HEADER    EXAMPLE\nATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00  0.00           N\nCONECT    1    2\n";
    let mut common = CommonRecords::new();
    for line in Lines::new(text) {
        common.observe(&line.expect("line"));
    }
    assert_eq!(common.lines.len(), 2);
    assert_eq!(fixed::record(common.lines[0].text), "HEADER");
    assert_eq!(fixed::record(common.lines[1].text), "CONECT");
}

#[test]
fn many_ragged_models_are_materialised_in_deposition_order() {
    let mut text = String::new();
    for model in 1..=128 {
        let atom = if model == 128 { " O  " } else { " N  " };
        let written = writeln!(
            text,
            "MODEL     {model:4}\nATOM      1 {atom} GLY A   1      1.000   1.000   1.000  1.00  0.00           N\nENDMDL"
        );
        assert!(written.is_ok());
    }
    let input = InputBuffer::from_bytes(text.into_bytes());
    let (structure, _) = crate::read(&input, &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("read failed: {findings:?}"));
    let Some(models) = structure.ragged_models() else {
        panic!("the final identity change must select ragged storage");
    };
    assert_eq!(models.len(), 128);
    let numbers: Vec<_> = structure
        .data()
        .models()
        .filter_map(pdbiox_core::structure::ModelRef::number)
        .collect();
    assert_eq!(numbers.first(), Some(&1));
    assert_eq!(numbers.last(), Some(&128));
}

#[test]
fn a_non_coordinate_atom_change_selects_independent_annotations() {
    let text = "MODEL        1\nATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00 10.00           N\nENDMDL\nMODEL        2\nATOM      1  N   GLY A   1      2.000   2.000   2.000  0.50 20.00           N\nENDMDL\n";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (structure, _) = crate::read(&input, &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("read failed: {findings:?}"));
    let Some(models) = structure.ragged_models() else {
        panic!("per-atom annotation changes must be ragged");
    };
    let occupancies: Vec<_> = models
        .iter()
        .filter_map(|model| model.data().atoms().next().and_then(AtomRef::occupancy))
        .collect();
    assert_eq!(occupancies, [1.0, 0.5]);
}

#[test]
fn ragged_models_receive_shared_metadata_and_connectivity() {
    use crate::PdbHeadersExt as _;

    let text = "HEADER    TEST                                      01-JAN-00   9XYZ\nCRYST1   10.000   11.000   12.000  90.00  90.00  90.00 P 1\nMODEL        1\nATOM      1  N   GLY A   1      1.000   1.000   1.000  1.00 10.00           N\nATOM      2  C   GLY A   1      2.000   2.000   2.000  1.00 10.00           C\nENDMDL\nMODEL        2\nATOM      1  O   GLY A   1      1.000   1.000   1.000  1.00 10.00           O\nATOM      2  C   GLY A   1      2.000   2.000   2.000  1.00 10.00           C\nENDMDL\nCONECT    1    2\n";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let (structure, _) = crate::read(&input, &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("read failed: {findings:?}"));
    let Some(models) = structure.ragged_models() else {
        panic!("identity-changing models must be ragged");
    };
    assert!(models.iter().all(|model| model.data().bonds.len() == 1));
    assert!(models.iter().all(|model| model.data().cell.is_some()));
    assert!(structure.pdb_headers().is_some());
}

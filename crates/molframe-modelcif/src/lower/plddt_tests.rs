use super::lower_plddt_annotation;
use crate::{ModelCifError, lower};
use molframe_core::annotation::PLDDT_ANNOTATION;
use molframe_core::column::Presence;
use molframe_core::io::{InputBuffer, ReadOptions};

const FIXTURE: &str = include_str!("../../tests/fixtures/modelcif_plddt.cif");

fn fixture() -> (molframe_core::Structure, crate::ModelCif) {
    let input = InputBuffer::from_bytes(FIXTURE.as_bytes().to_vec());
    let (document, structure, findings) =
        molframe_cif::read_with_document(&input, &ReadOptions::new())
            .expect("fixture structure lowers");
    assert!(findings.is_empty());
    let (model, findings) = lower(&document).expect("fixture ModelCIF lowers");
    assert!(findings.is_empty());
    (structure, model)
}

#[test]
fn residue_confidence_is_repeated_on_atoms_with_explicit_absence() {
    let (structure, model) = fixture();
    let column = lower_plddt_annotation(&structure, &model)
        .expect("pLDDT lowers")
        .expect("fixture has pLDDT");
    assert_eq!(column.len(), 10);
    assert_eq!(column.get(0), Some((92.5, Presence::Present)));
    assert_eq!(column.get(2), Some((92.5, Presence::Present)));
    assert_eq!(column.get(3), Some((0.0, Presence::Unknown)));
    assert_eq!(column.get(5), Some((0.0, Presence::Unknown)));
    assert_eq!(column.get(6), Some((55.0, Presence::Present)));
    assert_eq!(column.get(8), Some((55.0, Presence::Present)));
    assert_eq!(column.get(9), Some((0.0, Presence::Inapplicable)));
    assert!(structure.annotations().get(PLDDT_ANNOTATION).is_none());
}

#[test]
fn duplicate_residue_confidence_is_rejected() {
    let duplicate = FIXTURE.replace("1 A 3 SER 1 55.0", "1 A 1 ALA 1 80.0\n1 A 3 SER 1 55.0");
    let input = InputBuffer::from_bytes(duplicate.into_bytes());
    let (document, structure, _) =
        molframe_cif::read_with_document(&input, &ReadOptions::new()).expect("fixture lowers");
    let (model, _) = lower(&document).expect("ModelCIF lowers");
    assert!(matches!(
        lower_plddt_annotation(&structure, &model),
        Err(ModelCifError::DuplicatePlddtResidue { sequence: 1, .. })
    ));
}

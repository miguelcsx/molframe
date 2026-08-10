use super::*;
use pdbiox_core::{BondOrder, InputBuffer, ReadOptions};

const TWO_RESIDUES: &str = concat!(
    "ATOM      1  N   GLY A   1      11.000  12.000  13.000  1.00 20.00           N  \n",
    "ATOM      2  CA  GLY A   1      12.000  12.000  13.000  1.00 21.00           C  \n",
    "ATOM      3  N   ALA A   2      13.000  12.000  13.000  1.00 22.00           N  \n",
    "CONECT    1    2\n",
    "CONECT    2    3\n",
    "END\n",
);

#[test]
fn projection_preserves_hierarchy_chemistry_and_coordinates() {
    let input = InputBuffer::from_bytes(TWO_RESIDUES.as_bytes().to_vec());
    let result = pdbiox_pdb::read(&input, &ReadOptions::default());
    assert!(result.is_ok());
    let structure = match result {
        Ok((structure, _diagnostics)) => structure,
        Err(diagnostics) => panic!("reader failed: {diagnostics:?}"),
    };
    let export = TopologyExport::from_model(
        &structure,
        ModelIndex::new(0),
        pdbiox_core::contract::Namespace::Label,
    );
    assert!(export.is_ok());
    let export = match export {
        Ok(export) => export,
        Err(error) => panic!("projection failed: {error}"),
    };

    assert_eq!(export.chains[0].id, "A");
    assert_eq!(export.chains[0].residues, 0..2);
    assert_eq!(export.residues[1].name, "ALA");
    assert_eq!(export.residues[1].number, Some(2));
    assert_eq!(export.atoms[0].atomic_number, 7);
    assert!((export.atoms[0].mass - 14.007).abs() < f64::EPSILON);
    for (observed, expected) in export.atoms[2].position.iter().zip([13.0, 12.0, 13.0]) {
        assert!((*observed - expected).abs() < f32::EPSILON);
    }
    assert_eq!(export.bonds.len(), 2);
    assert_eq!(export.bonds[0].order, BondOrder::Unknown);
}

#[test]
fn unavailable_model_is_explicit() {
    let input = InputBuffer::from_bytes(TWO_RESIDUES.as_bytes().to_vec());
    let result = pdbiox_pdb::read(&input, &ReadOptions::default());
    let structure = match result {
        Ok((structure, _diagnostics)) => structure,
        Err(diagnostics) => panic!("reader failed: {diagnostics:?}"),
    };
    assert_eq!(
        TopologyExport::from_model(
            &structure,
            ModelIndex::new(1),
            pdbiox_core::contract::Namespace::Label
        ),
        Err(TopologyExportError::ModelUnavailable { model: 1 })
    );
}

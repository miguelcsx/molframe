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
fn projection_is_columnar_and_preserves_hierarchy() {
    let structure = structure();
    let batch = TopologyBatch::from_model(
        &structure,
        ModelIndex::new(0),
        pdbiox_core::contract::Namespace::Label,
    );
    let batch = match batch {
        Ok(batch) => batch,
        Err(error) => panic!("projection failed: {error}"),
    };
    assert_eq!(batch.string(batch.chain_ids[0]), Some("A"));
    assert_eq!(batch.chain_residue_offsets, [0, 2]);
    assert_eq!(batch.string(batch.residue_names[1]), Some("ALA"));
    assert_eq!(batch.residue_numbers[1], 2);
    assert_eq!(batch.atomic_numbers[0], 7);
    assert!((batch.masses[0] - 14.007).abs() < f64::EPSILON);
    assert!((batch.position_x[2] - 13.0).abs() < f32::EPSILON);
    assert!((batch.position_y[2] - 12.0).abs() < f32::EPSILON);
    assert!((batch.position_z[2] - 13.0).abs() < f32::EPSILON);
    assert_eq!(batch.bond_orders, [BondOrder::Unknown, BondOrder::Unknown]);
}

#[test]
fn unavailable_model_is_explicit() {
    assert_eq!(
        TopologyBatch::from_model(
            &structure(),
            ModelIndex::new(1),
            pdbiox_core::contract::Namespace::Label,
        ),
        Err(TopologyBatchError::ModelUnavailable { model: 1 })
    );
}

#[test]
fn columnar_projection_round_trips() {
    let source = structure();
    let batch = match TopologyBatch::from_model(
        &source,
        ModelIndex::new(0),
        pdbiox_core::contract::Namespace::Auth,
    ) {
        Ok(batch) => batch,
        Err(error) => panic!("projection failed: {error}"),
    };
    let imported = match batch.into_structure() {
        Ok(structure) => structure,
        Err(error) => panic!("import failed: {error}"),
    };
    assert_eq!(imported.chain_count(), 1);
    assert_eq!(imported.residue_count(), 2);
    assert_eq!(imported.atom_count(), 3);
    assert_eq!(imported.data().bonds.len(), 2);
    for (observed, expected) in imported.positions()[2].iter().zip([13.0, 12.0, 13.0]) {
        assert!((*observed - expected).abs() < f32::EPSILON);
    }
}

#[test]
fn import_rejects_a_misaligned_column() {
    let mut batch = match TopologyBatch::from_model(
        &structure(),
        ModelIndex::new(0),
        pdbiox_core::contract::Namespace::Label,
    ) {
        Ok(batch) => batch,
        Err(error) => panic!("projection failed: {error}"),
    };
    batch.position_z.pop();
    assert!(matches!(
        batch.into_structure(),
        Err(crate::TopologyImportError::ColumnLength {
            column: "position_z",
            ..
        })
    ));
}

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(TWO_RESIDUES.as_bytes().to_vec());
    match pdbiox_pdb::read(&input, &ReadOptions::default()) {
        Ok((structure, _diagnostics)) => structure,
        Err(diagnostics) => panic!("reader failed: {diagnostics:?}"),
    }
}

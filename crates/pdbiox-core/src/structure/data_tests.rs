use super::*;
use crate::coords::CoordinateBlock;

#[test]
fn an_empty_structure_reports_one_model_and_no_atoms() {
    let structure = Structure::new(StructureData::empty());
    assert_eq!(structure.atom_count(), 0);
    assert_eq!(structure.model_count(), 1);
    assert_eq!(structure.chain_count(), 0);
    assert!(structure.positions().is_empty());
}

#[test]
fn cloning_a_structure_shares_it_rather_than_copying_it() {
    let structure = Structure::new(StructureData::empty());
    let other = structure.clone();
    assert!(std::ptr::eq(structure.data(), other.data()));
}

#[test]
fn dense_models_report_one_block_each_and_ragged_ones_report_none() {
    let frames: Vec<CoordinateBlock> = (0..3).map(|_| (0..4).map(|_| [0.0; 3]).collect()).collect();
    let dense = CoordinateStore::Dense { frames };
    assert_eq!(dense.model_count(), 3);
    assert!(dense.is_dense());
    assert_eq!(
        dense.block(ModelIndex::new(2)).map(CoordinateBlock::len),
        Some(4)
    );
    assert!(dense.block(ModelIndex::new(3)).is_none());

    let ragged = CoordinateStore::Ragged {
        models: vec![Structure::new(StructureData::empty()); 2],
    };
    assert_eq!(ragged.model_count(), 2);
    assert!(!ragged.is_dense());
    assert!(ragged.block(ModelIndex::new(0)).is_none());
}

#[test]
fn a_unit_cube_with_right_angles_is_recognised_as_no_cell_at_all() {
    let placeholder = UnitCell {
        lengths: [1.0; 3],
        angles: [90.0; 3],
    };
    let real = UnitCell {
        lengths: [61.5, 61.5, 97.3],
        angles: [90.0, 90.0, 120.0],
    };
    assert!(placeholder.is_placeholder());
    assert!(!real.is_placeholder());
}

#[test]
fn positions_of_the_first_model_are_the_buffer_itself() {
    let mut data = StructureData::empty();
    data.coords = CoordinateStore::Single((0..3).map(|i| [i as f32, 0.0, 0.0]).collect());
    let structure = Structure::new(data);

    let first = structure.positions();
    assert_eq!(first.len(), 3);
    assert!(std::ptr::eq(first.as_ptr(), structure.positions().as_ptr()));
}

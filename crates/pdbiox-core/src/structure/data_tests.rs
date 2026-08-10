use super::*;
use crate::coords::CoordinateBlock;
use num_traits::ToPrimitive;

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
fn a_data_snapshot_clone_shares_atom_chunks() {
    let original = crate::structure::fixture::sample();
    let cloned = original.data().clone();
    assert!(Arc::ptr_eq(&original.data().chunks, &cloned.chunks));
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
    assert_eq!(ragged.ragged_models().map(<[Structure]>::len), Some(2));
}

#[test]
fn a_ragged_model_snapshot_is_its_own_single_model_structure() {
    let child = crate::structure::fixture::sample();
    let mut data = StructureData::empty();
    data.coords = CoordinateStore::Ragged {
        models: vec![child.clone()],
    };
    let ensemble = Structure::new(data);

    let Some((snapshot, local)) = ensemble.model_snapshot(ModelIndex::new(0)) else {
        panic!("ragged model missing")
    };
    assert!(std::ptr::eq(snapshot.data(), child.data()));
    assert_eq!(local, ModelIndex::new(0));
    assert!(ensemble.model_snapshot(ModelIndex::new(1)).is_none());
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
    data.coords = CoordinateStore::Single(
        (0..3)
            .map(|i| [i.to_f32().expect("small coordinate"), 0.0, 0.0])
            .collect(),
    );
    let structure = Structure::new(data);

    let first = structure.positions();
    assert_eq!(first.len(), 3);
    assert!(std::ptr::eq(first.as_ptr(), structure.positions().as_ptr()));
}

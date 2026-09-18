use super::{Unwrap, Wrap};
use crate::{FrameTransform, Timestep};
use molframe_core::structure::UnitCell;

fn timestep(positions: Vec<[f32; 3]>) -> Timestep {
    Timestep {
        positions,
        cell: Some(UnitCell {
            lengths: [10.0; 3],
            angles: [90.0; 3],
        }),
        ..Timestep::default()
    }
}

#[test]
fn bonded_molecule_is_made_whole_across_the_boundary() {
    let mut frame = timestep(vec![[9.8, 0.0, 0.0], [0.2, 0.0, 0.0], [0.6, 0.0, 0.0]]);
    Unwrap::molecules([(0, 1), (1, 2)])
        .apply(&mut frame)
        .unwrap_or_else(|error| panic!("unwrap failed: {error}"));
    assert!((frame.positions[1][0] - 10.2).abs() < 1.0e-5);
    assert!((frame.positions[2][0] - 10.6).abs() < 1.0e-5);
}

#[test]
fn atoms_and_groups_have_distinct_wrap_semantics() {
    let positions = vec![[9.8, 0.0, 0.0], [10.2, 0.0, 0.0]];
    let mut atoms = timestep(positions.clone());
    Wrap::atoms()
        .apply(&mut atoms)
        .unwrap_or_else(|error| panic!("atom wrap failed: {error}"));
    assert!(atoms.positions[1][0] < 1.0);

    let mut group = timestep(positions);
    Wrap::groups([[0, 1]])
        .apply(&mut group)
        .unwrap_or_else(|error| panic!("group wrap failed: {error}"));
    assert!((group.positions[1][0] - 0.2).abs() < 1.0e-5);
    assert!((group.positions[0][0] + 0.2).abs() < 1.0e-5);
}

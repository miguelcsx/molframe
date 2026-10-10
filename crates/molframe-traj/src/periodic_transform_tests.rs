use super::{TemporalUnwrap, Unwrap, Wrap};
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

#[test]
fn make_molecules_whole_is_an_alias_of_molecules() {
    let mut frame = timestep(vec![[9.8, 0.0, 0.0], [0.2, 0.0, 0.0]]);
    Unwrap::make_molecules_whole([(0, 1)])
        .apply(&mut frame)
        .unwrap_or_else(|error| panic!("unwrap failed: {error}"));
    assert!((frame.positions[1][0] - 10.2).abs() < 1.0e-5);
}

#[test]
fn a_molecule_crossing_the_boundary_has_a_continuous_centre() {
    // A two-atom molecule (1 A apart) drifting +0.4 A per frame along x in a
    // 10 A box; deposited positions are wrapped into [0, 10).
    let mut temporal = TemporalUnwrap::molecules([(0, 1)]);
    let mut whole_only = Unwrap::molecules([(0, 1)]);
    let mut previous_continuous: Option<f32> = None;
    let mut previous_whole: Option<f32> = None;
    let mut whole_jumped = false;
    for step in 0_u8..20 {
        let centre = 8.0 + 0.4 * f32::from(step);
        let wrap = |x: f32| x.rem_euclid(10.0);
        let positions = vec![
            [wrap(centre - 0.5), 0.0, 0.0],
            [wrap(centre + 0.5), 0.0, 0.0],
        ];
        let mut frame = timestep(positions.clone());
        temporal
            .apply(&mut frame)
            .unwrap_or_else(|error| panic!("temporal unwrap failed: {error}"));
        let mean = frame.positions[0][0].midpoint(frame.positions[1][0]);
        assert!(
            (frame.positions[1][0] - frame.positions[0][0] - 1.0).abs() < 1.0e-4,
            "molecule must stay whole"
        );
        if let Some(before) = previous_continuous {
            assert!(
                (mean - before - 0.4).abs() < 1.0e-4,
                "step {step}: {mean} vs {before}"
            );
        }
        previous_continuous = Some(mean);

        let mut whole = timestep(positions);
        whole_only
            .apply(&mut whole)
            .unwrap_or_else(|error| panic!("unwrap failed: {error}"));
        let whole_mean = whole.positions[0][0].midpoint(whole.positions[1][0]);
        if previous_whole.is_some_and(|before| (whole_mean - before).abs() > 5.0) {
            whole_jumped = true;
        }
        previous_whole = Some(whole_mean);
    }
    // The final continuous centre has left the box: 8 + 0.4 * 19 = 15.6.
    assert!(previous_continuous.is_some_and(|mean| (mean - 15.6).abs() < 1.0e-3));
    assert!(
        whole_jumped,
        "make-whole alone is expected to jump across the box"
    );
}

#[test]
fn temporal_unwrap_of_atoms_tracks_each_atom() {
    let mut transform = TemporalUnwrap::atoms();
    let mut first = timestep(vec![[9.9, 0.0, 0.0]]);
    let mut second = timestep(vec![[0.1, 0.0, 0.0]]);
    transform
        .apply(&mut first)
        .unwrap_or_else(|e| panic!("{e}"));
    transform
        .apply(&mut second)
        .unwrap_or_else(|e| panic!("{e}"));
    assert!((second.positions[0][0] - 10.1).abs() < 1.0e-5);
}

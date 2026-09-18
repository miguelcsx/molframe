use super::{PiStackingOptions, StackingKind, pi_stacking};
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => crate::chemistry_test_support::aromatic(
            &structure,
            0..u32::try_from(source.matches("ATOM ").count()).unwrap_or(0),
        ),
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn options() -> PiStackingOptions {
    PiStackingOptions {
        maximum_centre_distance: 6.0,
        maximum_parallel_angle: 30.0,
        minimum_perpendicular_angle: 60.0,
        plane_fit: molframe_geom::EigenOptions::standard(),
    }
}

// Ring 1 lies in the z = 0 plane.
const RING_ONE: &str = "\
ATOM 1 C CG PHE A 1 0 0 0\n\
ATOM 2 C CD1 PHE A 1 1 0 0\n\
ATOM 3 C CD2 PHE A 1 0 1 0\n";

#[test]
fn two_parallel_rings_stacked_face_to_face_are_reported() {
    // Ring 2 is the same pattern lifted 4 Å along z: parallel, centres 4 Å apart.
    let ring_two = "\
ATOM 4 C CG PHE A 2 0 0 4\n\
ATOM 5 C CD1 PHE A 2 1 0 4\n\
ATOM 6 C CD2 PHE A 2 0 1 4\n";
    let stacks = pi_stacking(
        &structure(&format!("{HEADER}{RING_ONE}{ring_two}")),
        options(),
    )
    .expect("options are valid");
    assert_eq!(stacks.len(), 1);
    assert_eq!((stacks[0].first.get(), stacks[0].second.get()), (0, 1));
    assert_eq!(stacks[0].kind, StackingKind::Parallel);
    assert!(stacks[0].angle < 1e-6, "angle {}", stacks[0].angle);
}

#[test]
fn perpendicular_rings_are_t_shaped() {
    // Ring 2 lies in the y = 0 plane, so its normal is perpendicular to ring 1's.
    let ring_two = "\
ATOM 4 C CG PHE A 2 0 0 4\n\
ATOM 5 C CD1 PHE A 2 1 0 4\n\
ATOM 6 C CD2 PHE A 2 0 0 5\n";
    let stacks = pi_stacking(
        &structure(&format!("{HEADER}{RING_ONE}{ring_two}")),
        options(),
    )
    .expect("options are valid");
    assert_eq!(stacks.len(), 1);
    assert_eq!(stacks[0].kind, StackingKind::TShaped);
    assert!(
        (stacks[0].angle - 90.0).abs() < 1e-4,
        "angle {}",
        stacks[0].angle
    );
}

#[test]
fn distant_rings_do_not_stack() {
    let ring_two = "\
ATOM 4 C CG PHE A 2 0 0 20\n\
ATOM 5 C CD1 PHE A 2 1 0 20\n\
ATOM 6 C CD2 PHE A 2 0 1 20\n";
    let stacks = pi_stacking(
        &structure(&format!("{HEADER}{RING_ONE}{ring_two}")),
        options(),
    )
    .expect("options are valid");
    assert!(stacks.is_empty());
}

#[test]
fn numerical_plane_failure_is_returned() {
    let mut policy = options();
    policy.plane_fit.relative_tolerance = 0.0;
    let error = pi_stacking(&structure(&format!("{HEADER}{RING_ONE}")), policy)
        .expect_err("invalid numerical controls must be returned");
    assert!(matches!(
        error,
        super::PiStackingError::Geometry(molframe_geom::EigenError::InvalidOptions)
    ));
}

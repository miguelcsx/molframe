use super::{
    PiStackingGeometry, PiStackingOptions, StackingKind, pi_stacking, pi_stacking_with_geometry,
};
use crate::ring_test_support::{Residue, cycle_bonds, hexagon, structure};
use molframe_core::structure::Structure;

const X: [f64; 3] = [1.0, 0.0, 0.0];
const Y: [f64; 3] = [0.0, 1.0, 0.0];
const Z: [f64; 3] = [0.0, 0.0, 1.0];

fn options() -> PiStackingOptions {
    PiStackingOptions {
        maximum_centre_distance: 6.0,
        maximum_parallel_angle: 30.0,
        minimum_perpendicular_angle: 60.0,
        plane_fit: molframe_geom::EigenOptions::standard(),
    }
}

/// Two single-ring residues, the second spanned by `u` and `v` about `centre`.
fn pair(centre: [f64; 3], u: [f64; 3], v: [f64; 3]) -> Structure {
    let mut bonds = cycle_bonds(0, 6);
    bonds.extend(cycle_bonds(6, 6));
    structure(
        &[
            Residue {
                name: "PHE",
                atoms: hexagon([0.0; 3], X, Y),
            },
            Residue {
                name: "PHE",
                atoms: hexagon(centre, u, v),
            },
        ],
        &bonds,
    )
}

fn stacks(structure: &Structure) -> Vec<super::PiStacking> {
    match pi_stacking(structure, options()) {
        Ok(table) => (0..table.len()).filter_map(|row| table.row(row)).collect(),
        Err(error) => panic!("valid options failed: {error}"),
    }
}

#[test]
fn rings_stacked_face_to_face_are_parallel() {
    let found = stacks(&pair([0.0, 0.0, 4.0], X, Y));
    assert_eq!(found.len(), 1);
    let stack = found[0];
    assert_eq!((stack.first.get(), stack.second.get()), (0, 1));
    assert_eq!(stack.kind, StackingKind::Parallel);
    assert!(stack.angle < 1e-6, "angle {}", stack.angle);
    assert!((stack.interplanar_separation - 4.0).abs() < 1e-4);
    assert!(stack.lateral_offset < 1e-3);
}

#[test]
fn displaced_stack_reports_the_offset() {
    let found = stacks(&pair([1.5, 0.0, 3.5], X, Y));
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].kind, StackingKind::Parallel);
    assert!((found[0].lateral_offset - 1.5).abs() < 1e-4);
    assert!((found[0].interplanar_separation - 3.5).abs() < 1e-4);
}

#[test]
fn coplanar_side_by_side_rings_are_not_parallel_stacks() {
    // Centroids 4.5 Å apart in the same plane: angle 0, distance within range.
    assert!(stacks(&pair([4.5, 0.0, 0.0], X, Y)).is_empty());
}

#[test]
fn parallel_planes_too_far_apart_or_too_offset_are_rejected() {
    // 5 Å separation is beyond the default 4.5 Å and beyond the 6 Å centre limit
    // once offset, and a 4 Å lateral shift at 3.2 Å separation exceeds 3.5 Å.
    assert!(stacks(&pair([4.0, 0.0, 3.2], X, Y)).is_empty());
    assert!(stacks(&pair([0.0, 0.0, 5.0], X, Y)).is_empty());
}

#[test]
fn edge_to_face_ring_pointing_at_the_face_is_t_shaped() {
    // Second ring stands in the x-z plane above the first ring's centre.
    let found = stacks(&pair([0.0, 0.0, 5.0], X, Z));
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].kind, StackingKind::TShaped);
    assert!((found[0].angle - 90.0).abs() < 1e-4);
}

#[test]
fn perpendicular_rings_beside_the_face_are_not_t_shaped() {
    // Perpendicular planes (normals z and x), but the centroid vector along y
    // lies in both planes: the rings meet edge to edge, not edge to face.
    assert!(stacks(&pair([0.0, 5.0, 0.0], Y, Z)).is_empty());
}

#[test]
fn distant_rings_do_not_stack() {
    assert!(stacks(&pair([0.0, 0.0, 20.0], X, Y)).is_empty());
}

#[test]
fn rings_of_one_residue_do_not_stack_on_each_other() {
    // A biaryl residue whose two rings are 4 Å apart face to face.
    let mut atoms = hexagon([0.0; 3], X, Y);
    atoms.extend(hexagon([0.0, 0.0, 4.0], X, Y));
    let mut bonds = cycle_bonds(0, 6);
    bonds.extend(cycle_bonds(6, 6));
    assert!(stacks(&structure(&[Residue { name: "BIP", atoms }], &bonds)).is_empty());
}

#[test]
fn each_ring_of_a_multi_ring_residue_is_paired_separately() {
    // Residue 1 holds two rings 4 Å above and 4 Å below residue 0's ring.
    let mut atoms = hexagon([0.0, 0.0, 4.0], X, Y);
    atoms.extend(hexagon([0.0, 0.0, -4.0], X, Y));
    let mut bonds = cycle_bonds(0, 6);
    bonds.extend(cycle_bonds(6, 6));
    bonds.extend(cycle_bonds(12, 6));
    let source = structure(
        &[
            Residue {
                name: "PHE",
                atoms: hexagon([0.0; 3], X, Y),
            },
            Residue { name: "BIP", atoms },
        ],
        &bonds.iter().map(|&(a, b)| (a, b)).collect::<Vec<_>>(),
    );
    let found = stacks(&source);
    assert_eq!(found.len(), 2);
    assert_ne!(found[0].second_atom, found[1].second_atom);
}

#[test]
fn periodic_image_ring_stacks_only_when_requested() {
    // Rings 14 Å apart along z in a 16 Å cell are 2 Å... use 18 Å cell: images 4 Å apart.
    let source = pair([0.0, 0.0, 14.0], X, Y);
    let mut data = source.data().clone();
    data.cell = Some(molframe_core::structure::UnitCell {
        lengths: [30.0, 30.0, 18.0],
        angles: [90.0; 3],
    });
    let source = Structure::new(data);
    let plain = pi_stacking(&source, options()).expect("valid options");
    assert!(plain.is_empty());
    let periodic = PiStackingGeometry {
        periodic: true,
        ..PiStackingGeometry::default()
    };
    let wrapped = pi_stacking_with_geometry(&source, options(), periodic).expect("valid cell");
    assert_eq!(wrapped.len(), 1);
    assert!((wrapped.row(0).expect("one stack").interplanar_separation - 4.0).abs() < 1e-3);
}

#[test]
fn periodic_geometry_rejects_missing_and_placeholder_cells() {
    let periodic = PiStackingGeometry {
        periodic: true,
        ..PiStackingGeometry::default()
    };
    let source = pair([0.0, 0.0, 4.0], X, Y);
    assert!(pi_stacking_with_geometry(&source, options(), periodic).is_err());
    let mut data = source.data().clone();
    data.cell = Some(molframe_core::structure::UnitCell {
        lengths: [1.0; 3],
        angles: [90.0; 3],
    });
    let source = Structure::new(data);
    assert!(pi_stacking_with_geometry(&source, options(), periodic).is_err());
}

#[test]
fn numerical_plane_failure_is_returned() {
    let mut policy = options();
    policy.plane_fit.relative_tolerance = 0.0;
    let error = pi_stacking(&pair([0.0, 0.0, 4.0], X, Y), policy)
        .expect_err("invalid numerical controls must be returned");
    assert!(matches!(
        error,
        super::PiStackingError::Geometry(molframe_geom::EigenError::InvalidOptions)
    ));
}

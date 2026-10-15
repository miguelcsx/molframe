use super::{CationPiOptions, cation_pi, cation_pi_periodic};
use crate::ring_test_support::{Residue, cycle_bonds, hexagon, structure};
use molframe_core::structure::Structure;

const X: [f64; 3] = [1.0, 0.0, 0.0];
const Y: [f64; 3] = [0.0, 1.0, 0.0];

/// A phenylalanine ring in the z = 0 plane centred on the origin, then a
/// lysine whose single charged atom sits at `cation`.
fn annotated_structure(cation: [f64; 3]) -> Structure {
    let source = structure(
        &[
            Residue {
                name: "PHE",
                atoms: hexagon([0.0; 3], X, Y),
            },
            Residue {
                name: "LYS",
                atoms: vec![cation],
            },
        ],
        &cycle_bonds(0, 6),
    );
    crate::chemistry_test_support::charges(&source, &[(6, 1)])
}

fn options() -> CationPiOptions {
    CationPiOptions {
        maximum_distance: 6.0,
        maximum_face_angle: 40.0,
        plane_fit: molframe_geom::EigenOptions::standard(),
    }
}

#[test]
fn a_lysine_amine_over_the_ring_face_is_a_cation_pi() {
    // NZ sits 4 Å above the ring centre, along the ring normal.
    let hits = cation_pi(&annotated_structure([0.0, 0.0, 4.0]), options()).expect("valid options");
    assert_eq!(hits.len(), 1);
    let hit = hits.row(0).expect("one cation-pi interaction");
    assert_eq!(hit.cation_residue.get(), 1);
    assert_eq!(hit.ring_residue.get(), 0);
    assert_eq!(hit.ring_atom.get(), 0);
    assert!((hit.distance - 4.0).abs() < 1e-4);
}

#[test]
fn a_cation_beside_the_ring_edge_is_not_counted() {
    let hits = cation_pi(&annotated_structure([5.0, 0.0, 0.0]), options()).expect("valid options");
    assert!(hits.is_empty());
}

#[test]
fn a_distant_cation_is_not_counted() {
    let hits = cation_pi(&annotated_structure([0.0, 0.0, 20.0]), options()).expect("valid options");
    assert!(hits.is_empty());
}

#[test]
fn a_cation_over_one_ring_of_a_fused_residue_names_that_ring() {
    // Residue 0 has two rings side by side in z = 0; the cation sits over the second.
    let mut atoms = hexagon([0.0; 3], X, Y);
    atoms.extend(hexagon([10.0, 0.0, 0.0], X, Y));
    let mut bonds = cycle_bonds(0, 6);
    bonds.extend(cycle_bonds(6, 6));
    let source = structure(
        &[
            Residue { name: "BIP", atoms },
            Residue {
                name: "LYS",
                atoms: vec![[10.0, 0.0, 4.0]],
            },
        ],
        &bonds,
    );
    let source = crate::chemistry_test_support::charges(&source, &[(12, 1)]);
    let hits = cation_pi(&source, options()).expect("valid options");
    assert_eq!(hits.len(), 1);
    assert_eq!(hits.row(0).expect("one hit").ring_atom.get(), 6);
}

#[test]
fn periodic_image_cation_is_found_only_when_requested() {
    // Direct distance 16 Å; through the 18 Å cell the image is 4 Å away.
    let source = annotated_structure([0.0, 0.0, 16.0]);
    let mut data = source.data().clone();
    data.cell = Some(molframe_core::structure::UnitCell {
        lengths: [30.0, 30.0, 20.0],
        angles: [90.0; 3],
    });
    let source = Structure::new(data);
    assert!(
        cation_pi(&source, options())
            .expect("valid options")
            .is_empty()
    );
    let hits = cation_pi_periodic(&source, options(), true).expect("valid cell");
    assert_eq!(hits.len(), 1);
    assert!((hits.row(0).expect("one hit").distance - 4.0).abs() < 1e-3);
}

#[test]
fn periodic_geometry_rejects_missing_and_placeholder_cells() {
    let source = annotated_structure([0.0, 0.0, 4.0]);
    assert!(cation_pi_periodic(&source, options(), true).is_err());
    let mut data = source.data().clone();
    data.cell = Some(molframe_core::structure::UnitCell {
        lengths: [1.0; 3],
        angles: [90.0; 3],
    });
    assert!(cation_pi_periodic(&Structure::new(data), options(), true).is_err());
}

#[test]
fn ring_fit_failure_is_not_hidden() {
    let mut policy = options();
    policy.plane_fit.maximum_sweeps = 0;
    let error = cation_pi(&annotated_structure([0.0, 0.0, 4.0]), policy)
        .expect_err("invalid ring-fit controls must be returned");
    assert!(matches!(
        error,
        super::CationPiError::RingGeometry(super::PiStackingError::Geometry(
            molframe_geom::EigenError::InvalidOptions
        ))
    ));
}

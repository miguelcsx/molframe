use super::{CationPiOptions, cation_pi};
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

// A phenylalanine ring in the z = 0 plane, centred near (0.33, 0.33, 0).
const RING: &str = "\
ATOM 1 C CG PHE A 1 0 0 0\n\
ATOM 2 C CD1 PHE A 1 1 0 0\n\
ATOM 3 C CD2 PHE A 1 0 1 0\n";

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn annotated_structure(cation: &str) -> Structure {
    let structure = crate::chemistry_test_support::aromatic(
        &structure(&format!("{HEADER}{RING}{cation}")),
        0..3,
    );
    crate::chemistry_test_support::charges(&structure, &[(3, 1)])
}

fn options() -> CationPiOptions {
    CationPiOptions {
        maximum_distance: 6.0,
        maximum_face_angle: 40.0,
        plane_fit: pdbiox_geom::EigenOptions::standard(),
    }
}

#[test]
fn a_lysine_amine_over_the_ring_face_is_a_cation_pi() {
    // NZ sits 4 Å above the ring centre, along the ring normal.
    let lysine = "ATOM 4 N NZ LYS A 2 0.33 0.33 4\n";
    let structure = annotated_structure(lysine);
    let hits = cation_pi(&structure, options()).expect("options are valid");
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].cation_residue.get(), 1);
    assert_eq!(hits[0].ring_residue.get(), 0);
    assert!((hits[0].distance - 4.0).abs() < 1e-4);
}

#[test]
fn a_cation_beside_the_ring_edge_is_not_counted() {
    // NZ level with the ring, off to the side: over the edge, not the face.
    let lysine = "ATOM 4 N NZ LYS A 2 5 0.33 0\n";
    let hits = cation_pi(&annotated_structure(lysine), options()).expect("options are valid");
    assert!(hits.is_empty());
}

#[test]
fn a_distant_cation_is_not_counted() {
    let lysine = "ATOM 4 N NZ LYS A 2 0.33 0.33 20\n";
    let hits = cation_pi(&annotated_structure(lysine), options()).expect("options are valid");
    assert!(hits.is_empty());
}

#[test]
fn ring_fit_failure_is_not_hidden() {
    let mut policy = options();
    policy.plane_fit.maximum_sweeps = 0;
    let error = cation_pi(
        &annotated_structure("ATOM 4 N NZ LYS A 2 0.33 0.33 4\n"),
        policy,
    )
    .expect_err("invalid ring-fit controls must be returned");
    assert!(matches!(
        error,
        super::CationPiError::RingGeometry(super::PiStackingError::Geometry(
            pdbiox_geom::EigenError::InvalidOptions
        ))
    ));
}

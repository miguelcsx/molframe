use super::aromatic_rings;
use crate::ring_test_support::{Residue, cycle_bonds, hexagon, structure};

const X: [f64; 3] = [1.0, 0.0, 0.0];
const Y: [f64; 3] = [0.0, 1.0, 0.0];

fn rings_of(structure: &molframe_core::structure::Structure) -> Vec<super::AromaticRing> {
    match aromatic_rings(structure, molframe_geom::EigenOptions::standard()) {
        Ok(rings) => rings,
        Err(error) => panic!("ring extraction failed: {error:?}"),
    }
}

#[test]
fn tryptophan_indole_yields_a_five_and_a_six_membered_ring() {
    // Atoms 0..6 form the benzene ring; 6, 7, 8 close the pyrrole ring on bond 0-1.
    let mut atoms = hexagon([0.0; 3], X, Y);
    atoms.extend([[2.3, -2.0, 0.0], [1.2, -3.0, 0.0], [0.0, -2.2, 0.0]]);
    let mut bonds = cycle_bonds(0, 6);
    bonds.extend([(1, 6), (6, 7), (7, 8), (8, 0)]);
    let trp = structure(&[Residue { name: "TRP", atoms }], &bonds);
    let rings = rings_of(&trp);
    let sizes: Vec<_> = rings.iter().map(|ring| ring.atoms.len()).collect();
    assert_eq!(sizes.len(), 2, "indole has two rings, got {sizes:?}");
    assert!(sizes.contains(&5) && sizes.contains(&6), "sizes {sizes:?}");
    assert!(rings.iter().all(|ring| ring.residue.get() == 0));
}

#[test]
fn naphthalene_has_two_rings_not_the_ten_atom_envelope() {
    let mut atoms = hexagon([0.0; 3], X, Y);
    atoms.extend([
        [2.4, -1.2, 0.0],
        [3.6, -0.5, 0.0],
        [3.6, 0.9, 0.0],
        [2.4, 1.6, 0.0],
    ]);
    let mut bonds = cycle_bonds(0, 6);
    bonds.extend([(1, 6), (6, 7), (7, 8), (8, 9), (9, 0)]);
    let rings = rings_of(&structure(&[Residue { name: "NAP", atoms }], &bonds));
    assert_eq!(rings.len(), 2);
}

#[test]
fn biphenyl_ligand_yields_two_rings_with_their_own_planes() {
    let mut atoms = hexagon([0.0; 3], X, Y);
    // Second ring twisted 60 degrees about the inter-ring axis (x).
    let twisted = [0.0, 0.5, 3.0_f64.sqrt() / 2.0];
    atoms.extend(hexagon([4.3, 0.0, 0.0], X, twisted));
    let mut bonds = cycle_bonds(0, 6);
    bonds.extend(cycle_bonds(6, 6));
    bonds.push((0, 6));
    let rings = rings_of(&structure(&[Residue { name: "BIP", atoms }], &bonds));
    assert_eq!(rings.len(), 2);
    assert!(rings.iter().all(|ring| ring.atoms.len() == 6));
    let cosine: f64 = (0..3)
        .map(|k| rings[0].normal[k] * rings[1].normal[k])
        .sum();
    let angle = cosine.abs().acos().to_degrees();
    assert!((angle - 60.0).abs() < 1e-3, "inter-ring angle {angle}");
    assert!((rings[0].centroid[0]).abs() < 1e-6 && (rings[1].centroid[0] - 4.3).abs() < 1e-6);
}

#[test]
fn atoms_without_bonds_form_no_rings() {
    let atoms = hexagon([0.0; 3], X, Y);
    assert!(rings_of(&structure(&[Residue { name: "PHE", atoms }], &[])).is_empty());
}

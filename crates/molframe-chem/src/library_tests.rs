use super::{cluster, describe, fingerprint, tanimoto};
use crate::mol::{
    MolAtom, MolAtomMetadata, MolBond, MolBondMetadata, MolRecord, MolVersion, Molecule,
};
use molframe_core::Element;

/// A record from atomic numbers and `(a, b, order)` bonds, over heavy atoms only.
fn record(z: &[u8], bonds: &[(usize, usize, u8)]) -> MolRecord {
    MolRecord {
        name: "test".into(),
        program: "".into(),
        comment: "".into(),
        version: MolVersion::V2000,
        molecule: Molecule {
            atoms: z
                .iter()
                .map(|z| MolAtom {
                    element: Element::from_atomic_number(*z),
                    position: [0.0; 3],
                })
                .collect(),
            bonds: bonds
                .iter()
                .map(|&(first, second, order)| MolBond {
                    first,
                    second,
                    order,
                })
                .collect(),
        },
        atom_metadata: vec![MolAtomMetadata::default(); z.len()],
        bond_metadata: vec![MolBondMetadata::default(); bonds.len()],
        properties: Vec::new(),
    }
}

/// Aspirin, C9H8O4, Kekulé, no explicit hydrogens.
fn aspirin() -> MolRecord {
    record(
        &[6, 6, 8, 8, 6, 6, 6, 6, 6, 6, 6, 8, 8],
        &[
            (0, 1, 1),
            (1, 2, 2),
            (1, 3, 1),
            (3, 4, 1),
            (4, 5, 2),
            (5, 6, 1),
            (6, 7, 2),
            (7, 8, 1),
            (8, 9, 2),
            (9, 4, 1),
            (9, 10, 1),
            (10, 11, 2),
            (10, 12, 1),
        ],
    )
}

/// Ibuprofen, C13H18O2.
fn ibuprofen() -> MolRecord {
    record(
        &[6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 8],
        &[
            (0, 1, 1),
            (2, 1, 1),
            (1, 3, 1),
            (3, 4, 1),
            (4, 5, 2),
            (5, 6, 1),
            (6, 7, 2),
            (7, 8, 1),
            (8, 9, 2),
            (9, 4, 1),
            (7, 10, 1),
            (10, 11, 1),
            (10, 12, 1),
            (12, 13, 2),
            (12, 14, 1),
        ],
    )
}

/// Caffeine, C8H10N4O2.
fn caffeine() -> MolRecord {
    // atoms: 0 N1, 1 C2, 2 N3, 3 C4, 4 C5, 5 C6, 6 N7, 7 C8, 8 N9, 9 O2, 10 O6, 11-13 methyls
    record(
        &[7, 6, 7, 6, 6, 6, 7, 6, 7, 8, 8, 6, 6, 6],
        &[
            (0, 1, 1),
            (1, 2, 1),
            (2, 3, 1),
            (3, 4, 2),
            (4, 5, 1),
            (5, 0, 1),
            (4, 6, 1),
            (6, 7, 1),
            (7, 8, 2),
            (8, 3, 1),
            (1, 9, 2),
            (5, 10, 2),
            (0, 11, 1),
            (2, 12, 1),
            (6, 13, 1),
        ],
    )
}

fn close(a: f64, b: f64, tolerance: f64) -> bool {
    (a - b).abs() < tolerance
}

#[test]
fn aspirin_has_the_properties_a_chemist_knows_it_by() {
    let d = describe(&aspirin());
    assert_eq!(d.heavy_atoms, 13);
    assert!(
        close(d.molecular_weight.unwrap(), 180.159, 0.05),
        "{:?}",
        d.molecular_weight
    );
    assert_eq!(d.hydrogen_bond_donors, Some(1));
    assert_eq!(d.hydrogen_bond_acceptors, 4);
    assert_eq!(d.rotatable_bonds, 3);
    assert_eq!(d.rings, 1);
    assert_eq!(d.aromatic_rings, 1);
    assert_eq!(d.formal_charge, 0);
    assert!(close(d.fraction_csp3.unwrap(), 1.0 / 9.0, 1e-9));
}

#[test]
fn ibuprofen_has_the_properties_a_chemist_knows_it_by() {
    let d = describe(&ibuprofen());
    assert!(
        close(d.molecular_weight.unwrap(), 206.285, 0.05),
        "{:?}",
        d.molecular_weight
    );
    assert_eq!(d.hydrogen_bond_donors, Some(1));
    assert_eq!(d.hydrogen_bond_acceptors, 2);
    assert_eq!(d.rotatable_bonds, 4);
    assert_eq!(d.rings, 1);
    assert_eq!(d.aromatic_rings, 1);
}

#[test]
fn caffeine_has_two_rings_no_donors_and_no_rotatable_bonds() {
    let d = describe(&caffeine());
    assert!(
        close(d.molecular_weight.unwrap(), 194.19, 0.05),
        "{:?}",
        d.molecular_weight
    );
    assert_eq!(d.hydrogen_bond_donors, Some(0));
    assert_eq!(d.hydrogen_bond_acceptors, 6);
    assert_eq!(d.rotatable_bonds, 0);
    assert_eq!(d.rings, 2);
    // the conservative rule may not call the fused rings aromatic, but never over-counts
    assert!(d.aromatic_rings <= 2);
}

#[test]
fn the_rule_of_five_says_which_rules_it_could_check() {
    let r = describe(&aspirin()).rule_of_five();
    assert_eq!((r.violations, r.checked), (0, 3));
}

#[test]
fn aromatic_typed_bonds_without_hydrogens_get_no_mass_rather_than_a_guess() {
    // benzene-like ring written with aromatic bonds and no hydrogens
    let ring = record(
        &[6; 6],
        &[
            (0, 1, 4),
            (1, 2, 4),
            (2, 3, 4),
            (3, 4, 4),
            (4, 5, 4),
            (5, 0, 4),
        ],
    );
    let d = describe(&ring);
    assert_eq!(d.molecular_weight, None);
    assert_eq!(d.hydrogen_bond_donors, None);
    assert_eq!(d.rings, 1);
    assert_eq!(d.aromatic_rings, 1);
}

#[test]
fn fused_rings_are_counted_as_the_smallest_set() {
    // naphthalene skeleton: two six-membered rings sharing a bond
    let n = record(
        &[6; 10],
        &[
            (0, 1, 2),
            (1, 2, 1),
            (2, 3, 2),
            (3, 4, 1),
            (4, 5, 2),
            (5, 0, 1),
            (4, 6, 1),
            (6, 7, 2),
            (7, 8, 1),
            (8, 9, 2),
            (9, 3, 1),
        ],
    );
    assert_eq!(describe(&n).rings, 2);
}

#[test]
fn a_molecule_is_identical_to_itself_and_a_close_relative_is_closer_than_a_stranger() {
    let (a, i, c) = (
        fingerprint(&aspirin()),
        fingerprint(&ibuprofen()),
        fingerprint(&caffeine()),
    );
    assert!(close(tanimoto(&a, &a), 1.0, 1e-12));
    assert!(close(tanimoto(&a, &c), tanimoto(&c, &a), 1e-12));
    // aspirin methyl ester: acid OH replaced by OCH3
    let mut ester = aspirin();
    ester.molecule.atoms.push(MolAtom {
        element: Element::from_atomic_number(6),
        position: [0.0; 3],
    });
    ester.molecule.bonds.push(MolBond {
        first: 12,
        second: 13,
        order: 1,
    });
    ester.atom_metadata.push(MolAtomMetadata::default());
    ester.bond_metadata.push(MolBondMetadata::default());
    let e = fingerprint(&ester);
    assert!(tanimoto(&a, &e) > tanimoto(&a, &c));
    assert!(tanimoto(&a, &e) > tanimoto(&a, &i) * 0.9);
    assert!(a.count() > 10);
}

#[test]
fn clusters_put_relatives_together_and_every_molecule_in_exactly_one() {
    let mut ester = aspirin();
    ester.molecule.atoms.push(MolAtom {
        element: Element::from_atomic_number(6),
        position: [0.0; 3],
    });
    ester.molecule.bonds.push(MolBond {
        first: 12,
        second: 13,
        order: 1,
    });
    ester.atom_metadata.push(MolAtomMetadata::default());
    ester.bond_metadata.push(MolBondMetadata::default());
    let prints = [
        fingerprint(&aspirin()),
        fingerprint(&ester),
        fingerprint(&caffeine()),
    ];
    let clusters = cluster(&prints, 0.4);
    let mut all: Vec<usize> = clusters.iter().flatten().copied().collect();
    all.sort_unstable();
    assert_eq!(all, vec![0, 1, 2]);
    assert!(
        clusters.iter().any(|c| c.contains(&0) && c.contains(&1)),
        "{clusters:?}"
    );
    assert!(clusters.iter().any(|c| c == &vec![2]), "{clusters:?}");
}

#[test]
fn a_threshold_of_one_leaves_every_molecule_alone_unless_identical() {
    let prints = [fingerprint(&aspirin()), fingerprint(&caffeine())];
    assert_eq!(cluster(&prints, 1.0).len(), 2);
    assert_eq!(cluster(&[], 0.5), Vec::<Vec<usize>>::new());
}

use super::*;
use crate::prelude::{ReadOptions, read_bytes};

const IDENTITY: [f32; 16] = [
    1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
];

fn hetatm(serial: u32, name: &str, element: &str, x: f32) -> String {
    format!(
        "HETATM{serial:>5} {name:<4} LIG A   1    {x:>8.3}{y:>8.3}{z:>8.3}  1.00  0.00          {element:>2}\n",
        y = if serial == 2 { 1.0 } else { 0.0 },
        z = 0.0,
    )
}

/// A three-atom chain C-C-O (so the two carbons are not equivalent), atoms named as given.
fn chain(names: [&str; 3], xs: [f32; 3]) -> Structure {
    let mut text = String::new();
    text.push_str(&hetatm(1, names[0], "C", xs[0]));
    text.push_str(&hetatm(2, names[1], "C", xs[1]));
    text.push_str(&hetatm(3, names[2], "O", xs[2]));
    text.push_str("CONECT    1    2\nCONECT    2    1    3\nCONECT    3    2\nEND\n");
    read_bytes(text.into_bytes(), Some("lig.pdb"), &ReadOptions::new())
        .map_or_else(|e| unreachable!("fixture is valid: {e}"), |(s, _)| s)
}

#[test]
fn the_heavy_atom_graph_of_a_ligand_is_a_component_with_its_bonds() {
    let ligand = chain(["C1", "C2", "O3"], [0.0, 1.5, 3.0]);
    let heavy = heavy_atom_component(&ligand).expect("heavy atoms");
    assert_eq!(heavy.atoms, [0, 1, 2]);
    assert_eq!(heavy.component.atoms.len(), 3);
    assert_eq!(heavy.component.bonds.len(), 2);
}

#[test]
fn identical_poses_have_no_distance_and_a_moved_one_has_the_shift() {
    let a = chain(["C1", "C2", "O3"], [0.0, 1.5, 3.0]);
    let b = chain(["C1", "C2", "O3"], [1.0, 2.5, 4.0]);
    assert_eq!(pose_rmsd(&a, &a, 1000), Some(0.0));
    let shifted = pose_rmsd(&a, &b, 1000).expect("same ligand");
    assert!((shifted - 1.0).abs() < 1e-5, "{shifted}");
}

#[test]
fn a_different_graph_is_not_the_same_ligand() {
    let a = chain(["C1", "C2", "O3"], [0.0, 1.5, 3.0]);
    let short = read_bytes(
        format!("{}END\n", hetatm(1, "C1", "C", 0.0)).into_bytes(),
        Some("one.pdb"),
        &ReadOptions::new(),
    )
    .map(|(s, _)| s)
    .expect("one atom");
    assert_eq!(pose_rmsd(&a, &short, 1000), None);
}

#[test]
fn a_placed_model_is_paired_by_name_then_by_order() {
    let reference = chain(["C1", "C2", "O3"], [0.0, 1.5, 3.0]);
    let renamed = chain(["X1", "X2", "X3"], [0.0, 1.5, 3.0]);
    let (by_name, pairing) =
        placed_rmsd(&reference, &reference, &IDENTITY, 1000).expect("compared");
    assert_eq!((by_name, pairing), (0.0, LigandPairing::ByName));
    let (by_order, pairing) = placed_rmsd(&reference, &renamed, &IDENTITY, 1000).expect("compared");
    assert_eq!((by_order, pairing), (0.0, LigandPairing::ByOrder));
}

#[test]
fn a_placement_that_carries_the_model_back_removes_the_distance() {
    let reference = chain(["C1", "C2", "O3"], [0.0, 1.5, 3.0]);
    let moved = chain(["C1", "C2", "O3"], [2.0, 3.5, 5.0]);
    let mut back = IDENTITY;
    back[12] = -2.0;
    let (rmsd, _) = placed_rmsd(&reference, &moved, &back, 1000).expect("compared");
    assert!(rmsd < 1e-5, "{rmsd}");
}

#[test]
fn an_unpairable_model_says_why() {
    let reference = chain(["C1", "C2", "O3"], [0.0, 1.5, 3.0]);
    let one = read_bytes(
        format!("{}END\n", hetatm(1, "C1", "C", 0.0)).into_bytes(),
        Some("one.pdb"),
        &ReadOptions::new(),
    )
    .map(|(s, _)| s)
    .expect("one atom");
    assert_eq!(
        placed_rmsd(&reference, &one, &IDENTITY, 1000).unwrap_err(),
        LigandCompareError::AtomCount {
            reference: 3,
            model: 1
        }
    );
}

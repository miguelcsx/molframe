use super::*;

#[test]
fn termini_missing_atoms_and_chain_breaks_are_absent_not_zero() {
    let residues = [residue(0.0, true), residue(2.0, false), residue(4.0, false)];
    let torsions = backbone_torsions(&residues);
    assert_eq!(torsions.len(), 3);
    assert!(torsions[0].phi.is_none());
    assert!(torsions[0].psi.is_some());
    assert!(torsions[0].omega.is_some());
    assert!(torsions[1].phi.is_some());
    assert!(torsions[1].psi.is_none());
    assert!(torsions[2].phi.is_none());
}

#[test]
fn a_missing_backbone_atom_invalidates_only_dependent_torsions() {
    let first = residue(0.0, true);
    let mut second = residue(2.0, false);
    second.alpha_carbon = None;
    let torsions = backbone_torsions(&[first, second]);
    assert!(torsions[0].psi.is_some());
    assert!(torsions[0].omega.is_none());
    assert!(torsions[1].phi.is_none());
}

fn residue(offset: f32, connected_to_next: bool) -> BackboneResidue {
    BackboneResidue {
        nitrogen: Some([offset, 0.0, 0.0]),
        alpha_carbon: Some([offset + 0.5, 1.0, 0.0]),
        carbonyl_carbon: Some([offset + 1.0, 1.0, 1.0]),
        connected_to_next,
    }
}

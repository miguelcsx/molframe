use super::*;
use molframe_core::AtomAnnotation;

const PEPTIDE: &str = "\
ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  GLY A   1       0.500   1.000   0.000  1.00 10.00           C
ATOM      3  C   GLY A   1       1.000   1.000   1.000  1.00 10.00           C
ATOM      4  N   ALA A   2       2.000   0.000   0.000  1.00 10.00           N
ATOM      5  CA  ALA A   2       2.500   1.000   0.000  1.00 10.00           C
ATOM      6  C   ALA A   2       3.000   1.000   1.000  1.00 10.00           C
ATOM      7  N   SER A   3       4.000   0.000   0.000  1.00 10.00           N
ATOM      8  CA  SER A   3       4.500   1.000   0.000  1.00 10.00           C
ATOM      9  C   SER A   3       5.000   1.000   1.000  1.00 10.00           C
CONECT    3    4
CONECT    6    7
END
";

#[test]
fn structure_torsions_follow_explicit_connectivity_and_keep_residue_indices() {
    let structure = match crate::read_bytes(
        PEPTIDE.as_bytes().to_vec(),
        Some("peptide.pdb"),
        &crate::ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture read failed: {findings:?}"),
    };
    let mut data = structure.engine().data().clone();
    let roles = [
        PolymerAtomRole::PROTEIN_NITROGEN,
        PolymerAtomRole::PROTEIN_ALPHA_CARBON,
        PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
    ];
    let role_values = (0..9)
        .map(|index| roles[index % roles.len()].code())
        .collect();
    data.annotations.insert(
        molframe_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(
            molframe_core::AnnotationColumn::from_values(role_values).expect("small column"),
        ),
    );
    let structure = Structure::new(data);
    let Ok(torsions) = structure_backbone_torsions(&structure) else {
        panic!("explicit backbone roles should be valid");
    };

    assert_eq!(torsions.len(), 3);
    assert_eq!(torsions[0].residue, ResidueIndex::new(0));
    assert!(torsions[0].torsions.phi.is_none());
    assert!(torsions[0].torsions.psi.is_some());
    assert!(torsions[1].torsions.phi.is_some());
    assert!(torsions[2].torsions.psi.is_none());
}

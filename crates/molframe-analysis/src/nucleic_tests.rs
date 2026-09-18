use super::nucleic_torsions;
use molframe_chem::PolymerAtomRole;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_core::{AnnotationColumn, AtomAnnotation};

// One guanosine with backbone atoms placed so β = P–O5'–C5'–C4' is 180°. It is
// the only residue, so α, ε and ζ (which reach into neighbours) are absent.
const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 P P G A 1 -0.5 1 0\n\
ATOM 2 O O5' G A 1 0 0 0\n\
ATOM 3 C C5' G A 1 1.33 0 0\n\
ATOM 4 C C4' G A 1 1.83 -1 0\n\
ATOM 5 C C3' G A 1 3.0 -1 0\n\
ATOM 6 O O3' G A 1 3.5 -2 0\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => with_roles(
            &structure,
            &[
                PolymerAtomRole::NUCLEIC_PHOSPHATE,
                PolymerAtomRole::NUCLEIC_O5,
                PolymerAtomRole::NUCLEIC_C5,
                PolymerAtomRole::NUCLEIC_C4,
                PolymerAtomRole::NUCLEIC_C3,
                PolymerAtomRole::NUCLEIC_O3,
            ],
        ),
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

#[test]
fn within_residue_torsions_are_measured_and_terminal_ones_are_absent() {
    let Ok(records) = nucleic_torsions(&structure()) else {
        panic!("explicit nucleic roles should be valid");
    };
    assert_eq!(records.len(), 1);
    let record = records[0];
    let Some(beta) = record.beta else {
        panic!("β should be defined");
    };
    assert!((beta.abs() - 180.0).abs() < 1e-3, "β was {beta}");
    assert!(record.alpha.is_none(), "no previous residue, so no α");
    assert!(record.epsilon.is_none(), "no next residue, so no ε");
    assert!(record.zeta.is_none());
}

#[test]
fn a_residue_without_a_sugar_is_skipped() {
    let source = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA ALA A 1 0 0 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let structure = with_roles(&structure, &[PolymerAtomRole::UNKNOWN]);
    let Ok(records) = nucleic_torsions(&structure) else {
        panic!("explicit non-nucleic role should be valid");
    };
    assert!(records.is_empty());
}

fn with_roles(structure: &Structure, roles: &[PolymerAtomRole]) -> Structure {
    let mut data = structure.data().clone();
    data.annotations.insert(
        molframe_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(
            AnnotationColumn::from_values(roles.iter().map(|role| role.code()).collect())
                .expect("small annotation column"),
        ),
    );
    Structure::new(data)
}

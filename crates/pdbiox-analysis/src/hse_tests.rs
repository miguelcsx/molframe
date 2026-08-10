use super::half_sphere_exposure;
use pdbiox_chem::PolymerAtomRole;
use pdbiox_core::Presence;
use pdbiox_core::annotation::{AnnotationColumn, AtomAnnotation};
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;
use pdbiox_spatial::SpatialBackend;

// Residue 1's side chain points along +x. Residue 2 sits on the +x side (upper),
// residue 3 on the −x side (lower).
const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CB ALA A 1 1 0 0\n\
ATOM 3 C CA ALA A 2 5 0 0\n\
ATOM 4 C CB ALA A 2 6 0 0\n\
ATOM 5 C CA ALA A 3 -5 0 0\n\
ATOM 6 C CB ALA A 3 -6 0 0\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => with_roles(&structure),
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn with_roles(structure: &Structure) -> Structure {
    let entries: Vec<_> = structure
        .data()
        .atoms()
        .map(|atom| {
            let role = match atom.name() {
                Some("CA") => PolymerAtomRole::PROTEIN_ALPHA_CARBON,
                Some("CB") => PolymerAtomRole::PROTEIN_BETA_CARBON,
                _ => PolymerAtomRole::UNKNOWN,
            };
            (role.code(), Presence::Present)
        })
        .collect();
    let mut data = structure.data().clone();
    let _ = data.annotations.insert(
        pdbiox_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(AnnotationColumn::from_entries(entries).expect("small column")),
    );
    Structure::new(data)
}

#[test]
fn a_residue_counts_neighbours_on_each_side_of_its_side_chain() {
    let Ok(exposure) = half_sphere_exposure(&structure(), 13.0, SpatialBackend::BruteForce) else {
        panic!("valid");
    };
    assert_eq!(exposure.len(), 3);
    // Residue 1 (index 0): residue 2 is up (+x), residue 3 is down (−x).
    assert_eq!(exposure[0].residue.get(), 0);
    assert_eq!((exposure[0].upper, exposure[0].lower), (1, 1));
}

#[test]
fn residues_without_a_beta_carbon_are_skipped() {
    // Glycine has no CB, so it never appears in the result.
    let source = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let Ok(exposure) =
        half_sphere_exposure(&with_roles(&structure), 13.0, SpatialBackend::BruteForce)
    else {
        panic!("valid");
    };
    assert!(exposure.is_empty());
}

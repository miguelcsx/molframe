use super::nonplanar_aromatic_rings;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn aromatic_structure(source: &str) -> Structure {
    let structure = structure(source);
    let mut data = structure.data().clone();
    data.annotations.insert(
        pdbiox_core::AROMATIC_ATOM_ANNOTATION,
        pdbiox_core::AtomAnnotation::Boolean(
            pdbiox_core::AnnotationColumn::from_values(vec![true; structure.atom_count() as usize])
                .expect("small annotation column"),
        ),
    );
    Structure::new(data)
}

fn flags(structure: &Structure) -> Vec<super::PlanarityFlag> {
    match nonplanar_aromatic_rings(
        structure,
        super::PlanarityOptions {
            maximum_deviation: 0.05,
            plane_fit: pdbiox_geom::EigenOptions::standard(),
        },
    ) {
        Ok(flags) => flags,
        Err(error) => panic!("planarity fixture failed: {error}"),
    }
}

// A phenylalanine ring laid flat in the z = 0 plane.
const FLAT_RING: &str = "\
ATOM 1 C CG PHE A 1 0.0 0.0 0\n\
ATOM 2 C CD1 PHE A 1 1.4 0.0 0\n\
ATOM 3 C CD2 PHE A 1 0.7 1.2 0\n\
ATOM 4 C CE1 PHE A 1 2.1 1.2 0\n\
ATOM 5 C CE2 PHE A 1 1.4 2.4 0\n\
ATOM 6 C CZ PHE A 1 0.0 1.2 0\n";

#[test]
fn a_flat_aromatic_ring_is_not_flagged() {
    let flags = flags(&aromatic_structure(&format!("{HEADER}{FLAT_RING}")));
    assert!(flags.is_empty(), "a flat ring should pass: {flags:?}");
}

#[test]
fn a_puckered_ring_is_flagged() {
    // The same ring with CZ lifted 1 Å out of the plane.
    let puckered = "\
ATOM 1 C CG PHE A 1 0.0 0.0 0\n\
ATOM 2 C CD1 PHE A 1 1.4 0.0 0\n\
ATOM 3 C CD2 PHE A 1 0.7 1.2 0\n\
ATOM 4 C CE1 PHE A 1 2.1 1.2 0\n\
ATOM 5 C CE2 PHE A 1 1.4 2.4 0\n\
ATOM 6 C CZ PHE A 1 0.0 1.2 1.0\n";
    let flags = flags(&aromatic_structure(&format!("{HEADER}{puckered}")));
    assert_eq!(flags.len(), 1);
    assert_eq!(flags[0].residue.get(), 0);
    assert!(flags[0].deviation > 0.05);
}

#[test]
fn a_non_aromatic_residue_is_ignored() {
    let source = format!(
        "{HEADER}\
ATOM 1 C CA ALA A 1 0.0 0.0 0\n\
ATOM 2 C CB ALA A 1 1.5 0.0 0\n"
    );
    let flags = flags(&structure(&source));
    assert!(flags.is_empty());
}

#[test]
fn invalid_plane_fit_controls_are_returned() {
    let error = nonplanar_aromatic_rings(
        &aromatic_structure(&format!("{HEADER}{FLAT_RING}")),
        super::PlanarityOptions {
            maximum_deviation: 0.05,
            plane_fit: pdbiox_geom::EigenOptions {
                relative_tolerance: 0.0,
                maximum_sweeps: 1,
            },
        },
    )
    .expect_err("invalid plane-fit controls must be returned");
    assert!(matches!(
        error,
        super::PlanarityError::Geometry(pdbiox_geom::EigenError::InvalidOptions)
    ));
}

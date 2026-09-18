use super::{DockQOptions, dockq};
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(body: &str) -> Structure {
    let input = InputBuffer::from_bytes(format!("{HEADER}{body}").into_bytes());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

// Receptor chain A and ligand chain B, each three non-collinear atoms, close
// enough to form an interface.
const NATIVE: &str = "\
ATOM 1 C CA GLY A 1 0 0 0\n\
ATOM 2 C CA GLY A 2 2 0 0\n\
ATOM 3 C CA GLY A 3 0 2 0\n\
ATOM 4 C CA GLY B 1 1 1 3\n\
ATOM 5 C CA GLY B 2 3 1 3\n\
ATOM 6 C CA GLY B 3 1 3 3\n";

#[test]
fn a_model_identical_to_the_native_scores_one() {
    let native = structure(NATIVE);
    let model = structure(NATIVE);
    let Ok(result) = dockq(&model, &native, "A", "B", options()) else {
        panic!("valid");
    };
    assert!((result.fnat - 1.0).abs() < 1e-12);
    assert!(result.ligand_rmsd < 1e-5);
    assert!(result.interface_rmsd < 1e-5);
    assert!((result.score - 1.0).abs() < 1e-6, "score {}", result.score);
}

#[test]
fn a_grossly_displaced_ligand_scores_near_zero() {
    let native = structure(NATIVE);
    // The same receptor, but the ligand shoved 50 Å away.
    let model = structure(
        "ATOM 1 C CA GLY A 1 0 0 0\n\
ATOM 2 C CA GLY A 2 2 0 0\n\
ATOM 3 C CA GLY A 3 0 2 0\n\
ATOM 4 C CA GLY B 1 51 1 3\n\
ATOM 5 C CA GLY B 2 53 1 3\n\
ATOM 6 C CA GLY B 3 51 3 3\n",
    );
    let Ok(result) = dockq(&model, &native, "A", "B", options()) else {
        panic!("valid");
    };
    assert!(result.fnat.abs() < 1e-12, "fnat {}", result.fnat);
    assert!(result.score < 0.23, "score {}", result.score);
}

fn options() -> DockQOptions {
    DockQOptions {
        contact_distance: 5.0,
        ligand_scale: 8.5,
        interface_scale: 1.5,
    }
}

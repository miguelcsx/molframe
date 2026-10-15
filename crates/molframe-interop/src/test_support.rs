//! Shared fixtures for the unit tests of this crate.

use molframe_core::Structure;
use molframe_core::io::{InputBuffer, ReadOptions};

const TWO_ATOMS: &str = "data_t\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 1 2 3\n\
ATOM 2 N N ALA A 1 4 5 6\n";

/// Reads mmCIF text, panicking with the findings when the fixture is invalid.
pub(crate) fn read_cif(text: &str) -> Structure {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

/// Two alanine atoms at (1, 2, 3) and (4, 5, 6).
pub(crate) fn two_atoms() -> Structure {
    read_cif(TWO_ATOMS)
}

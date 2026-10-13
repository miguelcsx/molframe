use super::match_atoms_by_name;
use molframe_core::io::{AmbiguousResidueBoundaryPolicy, InputBuffer, ReadOptions};
use molframe_core::structure::Structure;

fn structure(rows: &str) -> Structure {
    let header = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n\
_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
_atom_site.label_asym_id\n";
    // Chains are delimited by `label_asym_id`; repeat each row's chain there.
    let mut doubled = String::new();
    for row in rows.lines() {
        let chain = row.split_whitespace().nth(5).unwrap_or(".");
        doubled.push_str(row);
        doubled.push(' ');
        doubled.push_str(chain);
        doubled.push('\n');
    }
    let rows = doubled;
    let input = InputBuffer::from_bytes(format!("{header}{rows}").into_bytes());
    // A repeated atom name inside one residue is a deliberate fixture, not an error.
    let options = ReadOptions::new()
        .ambiguous_residue_boundary_policy(AmbiguousResidueBoundaryPolicy::InferFromFileOrder);
    match molframe_cif::read(&input, &options) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

#[test]
fn atoms_pair_by_chain_residue_and_name_whatever_their_order() {
    let reference = structure(
        "ATOM 1 N N GLY A 1 0 0 0\nATOM 2 C CA GLY A 1 1 0 0\nATOM 3 C CA ALA A 2 2 0 0\n",
    );
    // Extra leading atom and a different order shift every index.
    let target = structure(
        "ATOM 1 C CA ALA B 9 9 0 0\nATOM 2 C CA ALA A 2 2 0 0\nATOM 3 C CA GLY A 1 1 0 0\nATOM 4 N N GLY A 1 0 0 0\n",
    );
    assert_eq!(
        match_atoms_by_name(&reference, &target),
        [Some(3), Some(2), Some(1)]
    );
}

#[test]
fn an_atom_without_a_counterpart_or_a_repeat_maps_to_none() {
    let reference = structure(
        "ATOM 1 C CA GLY A 1 0 0 0\nATOM 2 C CB GLY A 1 1 0 0\nATOM 3 C CA GLY A 1 2 0 0\n",
    );
    let target = structure("ATOM 1 C CA GLY A 1 5 0 0\n");
    assert_eq!(
        match_atoms_by_name(&reference, &target),
        [Some(0), None, None]
    );
}

#[test]
fn the_same_residue_number_in_another_chain_is_not_a_match() {
    let reference = structure("ATOM 1 C CA GLY A 1 0 0 0\n");
    let target = structure("ATOM 1 C CA GLY B 1 0 0 0\n");
    assert_eq!(match_atoms_by_name(&reference, &target), [None]);
}

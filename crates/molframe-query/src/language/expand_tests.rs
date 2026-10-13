use super::same_residue;
use molframe_core::index::AtomIndex;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;
use std::fmt::Write as _;

/// Three chains of 40 residues with three atoms each.
fn structure() -> Structure {
    let mut source = String::from(
        "data_e\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    let mut serial = 0;
    for chain in ["A", "B", "C"] {
        for residue in 1..=40 {
            for (atom, element) in [("N", "N"), ("CA", "C"), ("C", "C")] {
                serial += 1;
                let _ = writeln!(
                    source,
                    "ATOM {serial} {element} {atom} ALA {chain} {residue} {serial} 0 0"
                );
            }
        }
    }
    molframe_cif::read(
        &InputBuffer::from_bytes(source.into_bytes()),
        &ReadOptions::new(),
    )
    .map_or_else(
        |findings| panic!("fixture failed: {findings:?}"),
        |(structure, _)| structure,
    )
}

/// Residue expansion by definition: every universe atom sharing a residue with a selected one.
fn reference(structure: &Structure, universe: &AtomSelection, selected: &[u32]) -> Vec<u32> {
    let residue_of = |atom: u32| {
        structure
            .atom(AtomIndex::new(atom))
            .and_then(molframe_core::structure::AtomRef::residue)
            .map(molframe_core::structure::ResidueRef::index)
    };
    let wanted: Vec<_> = selected
        .iter()
        .filter_map(|&atom| residue_of(atom))
        .collect();
    (0..structure.atom_count())
        .filter(|&atom| universe.contains(atom))
        .filter(|&atom| residue_of(atom).is_some_and(|residue| wanted.contains(&residue)))
        .collect()
}

#[test]
fn expansion_matches_the_definition_for_scattered_selections_and_universes() {
    let structure = structure();
    let count = structure.atom_count();
    let everything = AtomSelection::All(structure.atom_count());
    let half = AtomSelection::from_sorted((0..count).filter(|atom| atom % 2 == 0).collect());
    for stride in [1, 2, 7, 50, 119] {
        let selected: Vec<u32> = (0..count).step_by(stride).collect();
        let selection = AtomSelection::from_sorted(selected.clone());
        for universe in [&everything, &half] {
            let got: Vec<u32> = same_residue(&structure, universe, &selection)
                .iter()
                .collect();
            assert_eq!(
                got,
                reference(&structure, universe, &selected),
                "stride {stride}"
            );
        }
    }
}

#[test]
fn empty_selection_expands_to_nothing() {
    let structure = structure();
    let everything = AtomSelection::All(structure.atom_count());
    assert!(same_residue(&structure, &everything, &AtomSelection::Empty).is_empty());
}

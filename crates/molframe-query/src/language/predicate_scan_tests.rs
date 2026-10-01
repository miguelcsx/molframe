//! The residue- and chunk-level scans must select exactly what a per-atom
//! walk of the structure would.

use crate::Query;
use molframe_bench::{Sample, structure};
use molframe_core::contract::AnalysisPolicy;
use molframe_core::structure::{AtomRef, ChainRef, ResidueRef, Structure};

fn selected(structure: &Structure, text: &str) -> Vec<u32> {
    let query = match Query::compile(text) {
        Ok(query) => query,
        Err(findings) => panic!("`{text}` failed to compile: {findings:?}"),
    };
    let policy = AnalysisPolicy::default();
    match query.evaluate(structure, &policy, &crate::Groups::new(), None) {
        Ok(evaluation) => evaluation.selection.iter().collect(),
        Err(findings) => panic!("`{text}` failed to evaluate: {findings:?}"),
    }
}

/// The atoms for which `accepts` holds, found by asking each one.
fn per_atom(
    structure: &Structure,
    accepts: impl Fn(AtomRef<'_>, ChainRef<'_>) -> bool,
) -> Vec<u32> {
    let mut found = Vec::new();
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            for atom in residue.atoms() {
                if accepts(atom, chain) {
                    found.push(atom.index().get());
                }
            }
        }
    }
    found
}

#[test]
fn a_chain_predicate_selects_the_atoms_of_those_chains() {
    for sample in [Sample::Medium, Sample::Ensemble] {
        let structure = structure(sample);
        for chain in structure.data().chains().filter_map(ChainRef::label) {
            let expected = per_atom(&structure, |_, owner| owner.label() == Some(chain));
            assert_eq!(
                selected(&structure, &format!("label_chain {chain}")),
                expected
            );
        }
    }
}

#[test]
fn an_author_chain_predicate_selects_the_atoms_of_those_chains() {
    let structure = structure(Sample::Medium);
    for chain in structure.data().chains().filter_map(ChainRef::auth_label) {
        let expected = per_atom(&structure, |_, owner| owner.auth_label() == Some(chain));
        assert_eq!(
            selected(&structure, &format!("auth_chain {chain}")),
            expected
        );
    }
}

#[test]
fn a_residue_name_predicate_selects_the_atoms_of_those_residues() {
    let structure = structure(Sample::Medium);
    let expected = per_atom(&structure, |atom, _| {
        atom.residue().and_then(ResidueRef::name) == Some("HEM")
    });
    assert!(!expected.is_empty(), "the fixture carries heme");
    assert_eq!(selected(&structure, "label_resname HEM"), expected);
}

#[test]
fn a_residue_number_range_selects_the_atoms_of_those_residues() {
    let structure = structure(Sample::Medium);
    let expected = per_atom(&structure, |atom, _| {
        atom.residue()
            .and_then(ResidueRef::auth_seq_id)
            .is_some_and(|number| (10..=20).contains(&number))
    });
    assert!(!expected.is_empty());
    assert_eq!(selected(&structure, "auth_resid 10:20"), expected);
}

#[test]
fn a_thresholded_atom_column_selects_exactly_the_atoms_that_pass() {
    let structure = structure(Sample::Medium);
    let high = per_atom(&structure, |atom, _| {
        atom.b_factor().is_some_and(|b| f64::from(b) > 30.0)
    });
    assert!(!high.is_empty() && high.len() < structure.atom_count() as usize);
    assert_eq!(selected(&structure, "bfactor > 30"), high);
    let partial = per_atom(&structure, |atom, _| {
        atom.occupancy().is_some_and(|value| f64::from(value) < 1.0)
    });
    assert_eq!(selected(&structure, "occupancy < 1"), partial);
}

#[test]
fn a_scan_inside_a_narrowed_universe_stays_inside_it() {
    let structure = structure(Sample::Medium);
    let first = per_atom(&structure, |_, chain| chain.label() == Some("A"));
    let expected: Vec<u32> = per_atom(&structure, |atom, _| {
        atom.b_factor().is_some_and(|b| f64::from(b) > 30.0)
    })
    .into_iter()
    .filter(|atom| first.contains(atom))
    .collect();
    assert_eq!(
        selected(&structure, "label_chain A and bfactor > 30"),
        expected
    );
    let complement = per_atom(&structure, |_, chain| chain.label() != Some("A"));
    assert_eq!(selected(&structure, "not label_chain A"), complement);
}

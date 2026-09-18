//! CCD component-graph integrity checks.

use crate::{ComponentAtom, ComponentBond};
use molframe_core::{Code, Diagnostic};
use std::collections::BTreeSet;

pub(super) fn validate_component_topology(
    atoms: &[ComponentAtom],
    bonds: &[ComponentBond],
) -> Result<(), Vec<Diagnostic>> {
    let mut names = BTreeSet::new();
    let mut findings = Vec::new();
    for atom in atoms {
        if !names.insert(atom.name.as_ref()) {
            findings.push(
                Diagnostic::new(Code::E2005)
                    .in_category("chem_comp_atom")
                    .with_context("atom_id", atom.name.to_string()),
            );
        }
    }
    let mut endpoints = BTreeSet::new();
    for bond in bonds {
        if !names.contains(bond.atom_a.as_ref()) || !names.contains(bond.atom_b.as_ref()) {
            findings.push(
                Diagnostic::new(Code::E2005)
                    .in_category("chem_comp_bond")
                    .with_context("atom_id_1", bond.atom_a.to_string())
                    .with_context("atom_id_2", bond.atom_b.to_string()),
            );
            continue;
        }
        let pair = ordered_pair(bond);
        if pair.0 == pair.1 || !endpoints.insert(pair) {
            findings.push(
                Diagnostic::new(Code::E2005)
                    .in_category("chem_comp_bond")
                    .with_context("atom_id_1", bond.atom_a.to_string())
                    .with_context("atom_id_2", bond.atom_b.to_string()),
            );
        }
    }
    if findings.is_empty() {
        Ok(())
    } else {
        Err(findings)
    }
}

fn ordered_pair(bond: &ComponentBond) -> (&str, &str) {
    if bond.atom_a.as_ref() <= bond.atom_b.as_ref() {
        (bond.atom_a.as_ref(), bond.atom_b.as_ref())
    } else {
        (bond.atom_b.as_ref(), bond.atom_a.as_ref())
    }
}

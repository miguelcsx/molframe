//! Per-unit residue annotation: matches observed atoms to dictionary chemistry.
//!
//! Each worker keeps a private component cache so the shared resolver is only
//! consulted once per component id per thread; outputs are merged later in
//! file order, which keeps the result independent of thread count.

use super::atom_chemistry::{AtomChemistry, HbondRoles, observed_stereo};
use super::component_bonds::add_component_bonds;
use super::component_index::ComponentIndex;
use super::resolver::{Resolved, Resolver};
use super::tally::Tally;
use crate::{Component, ComponentKind};
use molframe_core::bond::BondRecord;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::structure::{AtomRef, ResidueRef};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;

/// Everything one unit of residues contributes, applied later in file order.
#[derive(Default)]
pub(super) struct UnitOutput {
    pub(super) kinds: BTreeMap<ComponentKind, usize>,
    /// Nucleotide residues carrying a 2'-oxygen (ribo) and those without (deoxy).
    pub(super) sugars: SugarCounts,
    pub(super) chemistry: Vec<(usize, AtomChemistry)>,
    pub(super) bonds: Vec<BondRecord>,
    pub(super) inconsistent: Tally,
    pub(super) missing_atoms: Tally,
    pub(super) missing: BTreeSet<Box<str>>,
}

/// Ribo- and deoxyribonucleotide residue counts of a unit.
#[derive(Clone, Copy, Default)]
pub(super) struct SugarCounts {
    pub(super) ribo: usize,
    pub(super) deoxy: usize,
}

impl SugarCounts {
    pub(super) fn merge(&mut self, other: Self) {
        self.ribo += other.ribo;
        self.deoxy += other.deoxy;
    }
}

/// Per-thread caches and scratch buffers, reused across units.
#[derive(Default)]
pub(super) struct Worker {
    local: HashMap<Box<str>, Option<Arc<Resolved>>>,
}

impl Worker {
    pub(super) fn run(
        &mut self,
        resolver: &Resolver<'_>,
        residues: &[ResidueRef<'_>],
    ) -> Result<UnitOutput, Diagnostic> {
        let mut output = UnitOutput::default();
        let mut atoms: Vec<AtomRef<'_>> = Vec::new();
        for residue in residues {
            self.annotate_residue(resolver, *residue, &mut atoms, &mut output)?;
        }
        Ok(output)
    }

    fn resolve(
        &mut self,
        resolver: &Resolver<'_>,
        id: &str,
    ) -> Result<Option<Arc<Resolved>>, Diagnostic> {
        if let Some(hit) = self.local.get(id) {
            return Ok(hit.clone());
        }
        let found = resolver.lookup(id)?;
        self.local.insert(id.into(), found.clone());
        Ok(found)
    }

    fn annotate_residue<'a>(
        &mut self,
        resolver: &Resolver<'_>,
        residue: ResidueRef<'a>,
        atoms: &mut Vec<AtomRef<'a>>,
        output: &mut UnitOutput,
    ) -> Result<(), Diagnostic> {
        atoms.clear();
        atoms.extend(residue.atoms());
        // Nearly every residue is one component, so group without a map; only
        // a microheterogeneous residue pays for the ordered split.
        let first = atoms.iter().find_map(|atom| atom.component_name());
        let uniform = atoms
            .iter()
            .all(|atom| atom.component_name().is_none() || atom.component_name() == first);
        if uniform {
            let Some(id) = first else {
                return Ok(());
            };
            let named: Vec<AtomRef<'a>> = atoms
                .iter()
                .copied()
                .filter(|atom| atom.component_name().is_some())
                .collect();
            return self.annotate_variant(resolver, id, &named, residue, output);
        }
        let mut variants: BTreeMap<&str, Vec<AtomRef<'a>>> = BTreeMap::new();
        for atom in atoms.iter() {
            let Some(component_id) = atom.component_name() else {
                continue;
            };
            variants.entry(component_id).or_default().push(*atom);
        }
        for (component_id, group) in variants {
            self.annotate_variant(resolver, component_id, &group, residue, output)?;
        }
        Ok(())
    }

    fn annotate_variant(
        &mut self,
        resolver: &Resolver<'_>,
        component_id: &str,
        atoms: &[AtomRef<'_>],
        residue: ResidueRef<'_>,
        output: &mut UnitOutput,
    ) -> Result<(), Diagnostic> {
        let Some(entry) = self.resolve(resolver, component_id)? else {
            output.missing.insert(component_id.into());
            return Ok(());
        };
        let component = &*entry.component;
        *output.kinds.entry(component.kind).or_insert(0) += 1;
        if component.kind == ComponentKind::Nucleotide {
            // The 2'-hydroxyl is what makes a nucleotide ribo; a definition
            // without ring atoms at all says nothing either way.
            let has = |name: &str| component.atom(name).is_some();
            if has("O2'") || has("O2*") {
                output.sugars.ribo += 1;
            } else if has("C2'") || has("C2*") {
                output.sugars.deoxy += 1;
            }
        }
        if inconsistent(&entry.index, atoms) {
            output.inconsistent.record(component_id, residue.index(), 0);
            return Ok(());
        }
        report_missing_atoms(component, atoms, residue, &mut output.missing_atoms);
        annotate_atoms(&entry.index, atoms, &mut output.chemistry);
        add_component_bonds(component, atoms, &mut output.bonds);
        Ok(())
    }
}

fn annotate_atoms(
    index: &ComponentIndex,
    atoms: &[AtomRef<'_>],
    output: &mut Vec<(usize, AtomChemistry)>,
) {
    let component = index.component();
    for atom in atoms {
        let Some(name) = atom.name() else {
            continue;
        };
        let Some(expected) = index.atom(name) else {
            continue;
        };
        output.push((
            atom.index().as_usize(),
            AtomChemistry {
                element: expected.element,
                carbonyl: index.is_carbonyl(name),
                kind: component.kind,
                aromatic: expected.aromatic,
                charge: expected.charge,
                roles: HbondRoles {
                    donor: index.is_donor(name),
                    acceptor: index.is_acceptor(name),
                    weak_acceptor: index.is_weak_acceptor(name),
                },
                stereo: observed_stereo(component, atoms, *atom, expected.stereo),
            },
        ));
    }
}

fn inconsistent(index: &ComponentIndex, atoms: &[AtomRef<'_>]) -> bool {
    atoms.iter().any(|atom| {
        // Hydrogen naming varies between deposition and preparation tools.
        // Unmatched hydrogens remain unannotated; matched heavy atoms retain
        // their dictionary chemistry and only confirmed bonds are installed.
        !atom
            .element()
            .is_some_and(molframe_core::Element::is_hydrogen)
            && atom.name().is_some_and(|name| index.atom(name).is_none())
    })
}

fn report_missing_atoms(
    component: &Component,
    atoms: &[AtomRef<'_>],
    residue: ResidueRef<'_>,
    tally: &mut Tally,
) {
    let observed: BTreeSet<_> = atoms.iter().filter_map(|atom| atom.name()).collect();
    let absent = component
        .atoms
        .iter()
        .filter(|atom| !atom.leaving && !observed.contains(atom.name.as_ref()))
        .count();
    if absent > 0 {
        tally.record(&component.id, residue.index(), absent);
    }
}

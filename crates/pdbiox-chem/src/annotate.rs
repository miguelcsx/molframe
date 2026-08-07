//! Applying component chemistry to an immutable structure snapshot.

use crate::{Component, ComponentKind, ComponentProvider};
use pdbiox_core::BondOrder;
use pdbiox_core::bond::{BondProvenance, BondRecord, BondTableBuilder};
use pdbiox_core::contract::DictionaryVersion;
use pdbiox_core::diagnostic::{Code, Diagnostic, Diagnostics};
use pdbiox_core::index::ChainIndex;
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};
use pdbiox_core::topology::PolymerKind;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// A chemically annotated snapshot and what could not be resolved.
#[derive(Clone, Debug)]
pub struct ChemistryReport {
    /// New immutable snapshot with connectivity and chain classification.
    pub structure: Structure,
    /// Non-fatal coverage and consistency findings.
    pub findings: Vec<Diagnostic>,
    /// Exact component dictionary version used.
    pub dictionary_version: DictionaryVersion,
}

/// Adds CCD internal bonds and checked polymer links without replacing file
/// connectivity.
///
/// # Errors
///
/// Returns a provider diagnostic when component storage itself fails. Unknown
/// components are successful reduced-coverage results and appear in `findings`.
pub fn apply_component_chemistry(
    structure: &Structure,
    provider: &dyn ComponentProvider,
) -> Result<ChemistryReport, Diagnostic> {
    let mut data = structure.data().clone();
    let mut bonds = BondTableBuilder::new();
    for bond in data.bonds.iter() {
        bonds.push(bond);
    }
    let mut state = Annotator {
        provider,
        components: BTreeMap::new(),
        missing: BTreeSet::new(),
        bonds,
        findings: Diagnostics::new(),
    };
    for chain in structure.data().chains() {
        let mut kinds = BTreeSet::new();
        let residues: Vec<_> = chain.residues().collect();
        for residue in &residues {
            state.annotate_residue(*residue, &mut kinds)?;
        }
        classify_chain(&mut data, chain.index(), &kinds);
        link_polymer(
            &residues,
            &state.components,
            &mut state.bonds,
            &mut state.findings,
        );
    }
    for component in state.missing {
        state
            .findings
            .push(Diagnostic::new(Code::W3201).with_context("component", component));
    }
    data.bonds = state.bonds.finish();
    for finding in pdbiox_core::structure::validate(&data) {
        state.findings.push(finding);
    }
    Ok(ChemistryReport {
        structure: Structure::new(data),
        findings: state.findings.finish(),
        dictionary_version: provider.version().clone(),
    })
}

struct Annotator<'a> {
    provider: &'a dyn ComponentProvider,
    components: BTreeMap<Box<str>, Arc<Component>>,
    missing: BTreeSet<Box<str>>,
    bonds: BondTableBuilder,
    findings: Diagnostics,
}

impl Annotator<'_> {
    fn annotate_residue<'a>(
        &mut self,
        residue: ResidueRef<'a>,
        kinds: &mut BTreeSet<ComponentKind>,
    ) -> Result<(), Diagnostic> {
        let mut variants: BTreeMap<&str, Vec<AtomRef<'a>>> = BTreeMap::new();
        for atom in residue.atoms() {
            let Some(component_id) = atom.component_name() else {
                continue;
            };
            variants.entry(component_id).or_default().push(atom);
        }
        for (component_id, atoms) in variants {
            let component = match self.components.get(component_id) {
                Some(component) => Some(Arc::clone(component)),
                None => self.provider.get(component_id)?,
            };
            let Some(component) = component else {
                self.missing.insert(component_id.into());
                continue;
            };
            self.components
                .entry(component_id.into())
                .or_insert_with(|| Arc::clone(&component));
            kinds.insert(component.kind);
            if inconsistent(&component, &atoms) {
                self.findings.push(
                    Diagnostic::new(Code::W3202)
                        .with_context("component", component_id)
                        .with_context("residue", residue.index().to_string()),
                );
                continue;
            }
            report_missing_atoms(&component, &atoms, residue, &mut self.findings);
            add_component_bonds(&component, &atoms, &mut self.bonds);
        }
        Ok(())
    }
}

fn inconsistent(component: &Component, atoms: &[AtomRef<'_>]) -> bool {
    atoms.iter().any(|atom| {
        atom.name()
            .is_some_and(|name| component.atom(name).is_none())
    })
}

fn report_missing_atoms(
    component: &Component,
    atoms: &[AtomRef<'_>],
    residue: ResidueRef<'_>,
    findings: &mut Diagnostics,
) {
    let observed: BTreeSet<_> = atoms.iter().filter_map(|atom| atom.name()).collect();
    let absent = component
        .atoms
        .iter()
        .filter(|atom| !atom.leaving && !observed.contains(atom.name.as_ref()))
        .count();
    if absent > 0 {
        findings.push(
            Diagnostic::new(Code::W3301)
                .with_context("component", component.id.to_string())
                .with_context("residue", residue.index().to_string())
                .with_context("missing atoms", absent.to_string()),
        );
    }
}

fn add_component_bonds(
    component: &Component,
    atoms: &[AtomRef<'_>],
    output: &mut BondTableBuilder,
) {
    for bond in component.bonds.iter() {
        for atom_a in atoms
            .iter()
            .filter(|atom| atom.name() == Some(&bond.atom_a))
        {
            for atom_b in atoms
                .iter()
                .filter(|atom| atom.name() == Some(&bond.atom_b))
            {
                if alt_compatible(*atom_a, *atom_b) {
                    output.push(BondRecord {
                        atom_a: atom_a.index(),
                        atom_b: atom_b.index(),
                        order: bond.order,
                        provenance: BondProvenance::ChemicalComponentDictionary,
                    });
                }
            }
        }
    }
}

fn alt_compatible(atom_a: AtomRef<'_>, atom_b: AtomRef<'_>) -> bool {
    match (atom_a.alt_id(), atom_b.alt_id()) {
        (Some(a), Some(b)) => a.is_blank() || b.is_blank() || a == b,
        _ => true,
    }
}

fn classify_chain(
    data: &mut pdbiox_core::StructureData,
    chain: ChainIndex,
    kinds: &BTreeSet<ComponentKind>,
) {
    let kind = if kinds.iter().all(|kind| *kind == ComponentKind::AminoAcid) && !kinds.is_empty() {
        PolymerKind::Protein
    } else if kinds.iter().all(|kind| *kind == ComponentKind::Nucleotide) && !kinds.is_empty() {
        PolymerKind::NucleicHybrid
    } else if kinds.iter().all(|kind| *kind == ComponentKind::Saccharide) && !kinds.is_empty() {
        PolymerKind::Saccharide
    } else {
        return;
    };
    data.topology.chains.set_polymer_kind(chain, kind);
}

fn link_polymer(
    residues: &[ResidueRef<'_>],
    components: &BTreeMap<Box<str>, Arc<Component>>,
    bonds: &mut BondTableBuilder,
    findings: &mut Diagnostics,
) {
    for pair in residues.windows(2) {
        let [left, right] = pair else {
            continue;
        };
        let Some(left_component) = left.name().and_then(|name| components.get(name)) else {
            continue;
        };
        let Some(right_component) = right.name().and_then(|name| components.get(name)) else {
            continue;
        };
        let Some((left_name, right_name)) =
            linkage_names(left_component.kind, right_component.kind)
        else {
            continue;
        };
        let (Some(atom_a), Some(atom_b)) = (left.atom(left_name), right.atom(right_name)) else {
            continue;
        };
        if within_linkage(atom_a, atom_b) {
            bonds.push(BondRecord {
                atom_a: atom_a.index(),
                atom_b: atom_b.index(),
                order: BondOrder::Polymeric,
                provenance: BondProvenance::ChemicalComponentDictionary,
            });
        } else {
            findings.push(
                Diagnostic::new(Code::W3302)
                    .with_context("left residue", left.index().to_string())
                    .with_context("right residue", right.index().to_string()),
            );
        }
    }
}

fn linkage_names(
    left: ComponentKind,
    right: ComponentKind,
) -> Option<(&'static str, &'static str)> {
    match (left, right) {
        (ComponentKind::AminoAcid, ComponentKind::AminoAcid) => Some(("C", "N")),
        (ComponentKind::Nucleotide, ComponentKind::Nucleotide) => Some(("O3'", "P")),
        _ => None,
    }
}

fn within_linkage(atom_a: AtomRef<'_>, atom_b: AtomRef<'_>) -> bool {
    let (Some(a), Some(b)) = (atom_a.position(), atom_b.position()) else {
        return false;
    };
    let squared: f32 = a
        .iter()
        .zip(b)
        .map(|(left, right)| (left - right).powi(2))
        .sum();
    squared <= 2.1f32.powi(2)
}

#[cfg(test)]
#[path = "annotate_tests.rs"]
mod tests;

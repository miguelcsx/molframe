//! Applying component chemistry to an immutable structure snapshot.

mod atom_chemistry;

use crate::{Component, ComponentKind, ComponentProvider};
use atom_chemistry::{
    AtomChemistry, attach_annotations, is_hydrogen_bond_acceptor, is_hydrogen_bond_donor,
    observed_stereo,
};
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
    /// Explicit rule used to create inter-residue polymer bonds.
    pub polymer_link_policy: PolymerLinkPolicy,
}

/// Explicit policy for constructing bonds between consecutive polymer components.
#[derive(Clone, Debug, PartialEq)]
pub enum PolymerLinkPolicy {
    /// Do not construct inter-residue bonds.
    Disabled,
    /// Link CCD-declared atoms only when they are within this distance.
    Explicit {
        /// Inclusive distance in ångström.
        angstrom: f32,
        /// Ordered component/atom attachment rules.
        rules: Arc<[PolymerLinkRule]>,
    },
}

/// One explicit inter-component attachment rule.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolymerLinkRule {
    /// Component role on the earlier residue.
    pub left_kind: ComponentKind,
    /// Component-local attachment atom on the earlier residue.
    pub left_atom: Box<str>,
    /// Component role on the later residue.
    pub right_kind: ComponentKind,
    /// Component-local attachment atom on the later residue.
    pub right_atom: Box<str>,
}

impl PolymerLinkRule {
    /// Constructs a named attachment rule without interpreting atom names.
    #[must_use]
    pub fn new(
        left_kind: ComponentKind,
        left_atom: impl Into<Box<str>>,
        right_kind: ComponentKind,
        right_atom: impl Into<Box<str>>,
    ) -> Self {
        Self {
            left_kind,
            left_atom: left_atom.into(),
            right_kind,
            right_atom: right_atom.into(),
        }
    }
}

impl PolymerLinkPolicy {
    /// Constructs an explicit distance-and-attachment policy.
    #[must_use]
    pub fn explicit(angstrom: f32, rules: impl Into<Arc<[PolymerLinkRule]>>) -> Self {
        Self::Explicit {
            angstrom,
            rules: rules.into(),
        }
    }
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
    polymer_link_policy: PolymerLinkPolicy,
) -> Result<ChemistryReport, Diagnostic> {
    validate_link_policy(&polymer_link_policy)?;
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
        chemistry: vec![None; structure.atom_count() as usize],
    };
    for chain in structure.data().chains() {
        let mut kinds = BTreeSet::new();
        let residues: Vec<_> = chain.residues().collect();
        for residue in &residues {
            state.annotate_residue(*residue, &mut kinds)?;
        }
        classify_chain(&mut data, chain.index(), &kinds)?;
        link_polymer(
            &residues,
            &state.components,
            &mut state.bonds,
            &mut state.findings,
            &polymer_link_policy,
        );
    }
    for component in &state.missing {
        state
            .findings
            .push(Diagnostic::new(Code::W3201).with_context("component", component.to_string()));
    }
    attach_annotations(&mut data, &state.chemistry)?;
    data.bonds = state.bonds.finish();
    for finding in pdbiox_core::structure::validate(&data) {
        state.findings.push(finding);
    }
    Ok(ChemistryReport {
        structure: Structure::new(data),
        findings: state.findings.finish(),
        dictionary_version: provider.version().clone(),
        polymer_link_policy,
    })
}

fn validate_link_policy(policy: &PolymerLinkPolicy) -> Result<(), Diagnostic> {
    match policy {
        PolymerLinkPolicy::Disabled => Ok(()),
        PolymerLinkPolicy::Explicit { angstrom, rules }
            if angstrom.is_finite() && *angstrom > 0.0 && !rules.is_empty() =>
        {
            Ok(())
        }
        PolymerLinkPolicy::Explicit { .. } => Err(Diagnostic::new(Code::E4003).with_context(
            "required",
            "finite positive polymer-link distance and attachment rules",
        )),
    }
}

struct Annotator<'a> {
    provider: &'a dyn ComponentProvider,
    components: BTreeMap<Box<str>, Arc<Component>>,
    missing: BTreeSet<Box<str>>,
    bonds: BondTableBuilder,
    findings: Diagnostics,
    chemistry: Vec<Option<AtomChemistry>>,
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
            self.annotate_atoms(&component, residue, &atoms);
            add_component_bonds(&component, &atoms, &mut self.bonds);
        }
        Ok(())
    }

    fn annotate_atoms(
        &mut self,
        component: &Component,
        residue: ResidueRef<'_>,
        atoms: &[AtomRef<'_>],
    ) {
        for atom in atoms {
            let Some(name) = atom.name() else {
                continue;
            };
            let Some(expected) = component.atom(name) else {
                continue;
            };
            let Some(slot) = self.chemistry.get_mut(atom.index().as_usize()) else {
                continue;
            };
            *slot = Some(AtomChemistry {
                kind: component.kind,
                aromatic: expected.aromatic,
                charge: expected.charge,
                donor: is_hydrogen_bond_donor(component, name),
                acceptor: is_hydrogen_bond_acceptor(component, name),
                stereo: observed_stereo(component, residue, name, expected.stereo),
            });
        }
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
) -> Result<(), Diagnostic> {
    let kind = if kinds.iter().all(|kind| *kind == ComponentKind::AminoAcid) && !kinds.is_empty() {
        PolymerKind::Protein
    } else if kinds.iter().all(|kind| *kind == ComponentKind::Nucleotide) && !kinds.is_empty() {
        PolymerKind::NucleicHybrid
    } else if kinds.iter().all(|kind| *kind == ComponentKind::Saccharide) && !kinds.is_empty() {
        PolymerKind::Saccharide
    } else {
        return Ok(());
    };
    data.topology
        .chains
        .set_polymer_kind(chain, kind)
        .map_err(|error| Diagnostic::new(Code::E3001).with_context("cause", error.to_string()))
}

fn link_polymer(
    residues: &[ResidueRef<'_>],
    components: &BTreeMap<Box<str>, Arc<Component>>,
    bonds: &mut BondTableBuilder,
    findings: &mut Diagnostics,
    policy: &PolymerLinkPolicy,
) {
    let PolymerLinkPolicy::Explicit { angstrom, rules } = policy else {
        return;
    };
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
        let Some(rule) = rules.iter().find(|rule| {
            rule.left_kind == left_component.kind && rule.right_kind == right_component.kind
        }) else {
            continue;
        };
        let (Some(atom_a), Some(atom_b)) =
            (left.atom(&rule.left_atom), right.atom(&rule.right_atom))
        else {
            continue;
        };
        if within_linkage(atom_a, atom_b, *angstrom) {
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

fn within_linkage(atom_a: AtomRef<'_>, atom_b: AtomRef<'_>, max_distance: f32) -> bool {
    let (Some(a), Some(b)) = (atom_a.position(), atom_b.position()) else {
        return false;
    };
    let squared: f32 = a
        .iter()
        .zip(b)
        .map(|(left, right)| (left - right).powi(2))
        .sum();
    squared <= max_distance.powi(2)
}

#[cfg(test)]
#[path = "annotate_tests.rs"]
mod tests;

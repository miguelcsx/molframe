//! Applying component chemistry to an immutable structure snapshot.

mod atom_chemistry;
mod component_index;
mod linkage;

use component_index::ComponentIndex;
use linkage::link_polymer;

use crate::{Component, ComponentKind, ComponentProvider};
pub use atom_chemistry::HBOND_WEAK_ACCEPTOR_ANNOTATION;
use atom_chemistry::{AtomChemistry, attach_annotations, observed_stereo};
use molframe_core::bond::{BondProvenance, BondRecord, BondTableBuilder};
use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::{Code, Diagnostic, Diagnostics};
use molframe_core::index::ChainIndex;
use molframe_core::structure::{AtomRef, ResidueRef, Structure};
use molframe_core::topology::PolymerKind;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

/// Where the per-atom formal charges written by the pass come from.
///
/// A component dictionary records the charge of the free, isolated component in
/// its reference protonation state. It is not the charge the atom carries at the
/// pH or in the environment of the modelled structure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum FormalChargeSource {
    /// Copied from the component definition, with no protonation modelling.
    ComponentDefault,
}

/// Provenance for the native component-chemistry pass.
#[derive(Clone, Debug, PartialEq)]
pub struct ChemistryProvenance {
    /// Origin of the formal-charge annotation. Consumers that need the real
    /// protonation state must not treat a `ComponentDefault` charge as one.
    pub charge_source: FormalChargeSource,
    /// Stable identity of the native pipeline.
    pub algorithm: &'static str,
    /// Exact component dictionary version used for lookup.
    pub dictionary_version: DictionaryVersion,
    /// Explicit policy used for inter-residue polymer bonds.
    pub polymer_link_policy: PolymerLinkPolicy,
}

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
    /// Input provenance for the native chemistry pipeline.
    pub provenance: ChemistryProvenance,
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
        indices: BTreeMap::new(),
        missing: BTreeSet::new(),
        bonds,
        findings: Diagnostics::new(),
        chemistry: vec![None; structure.atom_count() as usize],
    };
    for chain in structure.data().chains() {
        let mut kinds = BTreeMap::new();
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
    for finding in molframe_core::structure::validate(&data) {
        state.findings.push(finding);
    }
    Ok(ChemistryReport {
        structure: Structure::new(data),
        findings: state.findings.finish(),
        dictionary_version: provider.version().clone(),
        polymer_link_policy: polymer_link_policy.clone(),
        provenance: ChemistryProvenance {
            charge_source: FormalChargeSource::ComponentDefault,
            algorithm: "component-chemistry",
            dictionary_version: provider.version().clone(),
            polymer_link_policy,
        },
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
    /// Name index and hydrogen-bonding roles, built once per component.
    ///
    /// Both are pure functions of the component, so every residue of the same
    /// component shares them instead of recomputing a quadratic bond walk.
    indices: BTreeMap<Box<str>, Arc<ComponentIndex>>,
    missing: BTreeSet<Box<str>>,
    bonds: BondTableBuilder,
    findings: Diagnostics,
    chemistry: Vec<Option<AtomChemistry>>,
}

impl Annotator<'_> {
    fn annotate_residue<'a>(
        &mut self,
        residue: ResidueRef<'a>,
        kinds: &mut BTreeMap<ComponentKind, usize>,
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
            let index = Arc::clone(
                self.indices
                    .entry(component_id.into())
                    .or_insert_with(|| Arc::new(ComponentIndex::build(Arc::clone(&component)))),
            );
            *kinds.entry(component.kind).or_insert(0) += 1;
            if inconsistent(&index, &atoms) {
                self.findings.push(
                    Diagnostic::new(Code::W3202)
                        .with_context("component", component_id)
                        .with_context("residue", residue.index().to_string()),
                );
                continue;
            }
            report_missing_atoms(&component, &atoms, residue, &mut self.findings);
            self.annotate_atoms(&index, &atoms);
            add_component_bonds(&component, &atoms, &mut self.bonds);
        }
        Ok(())
    }

    fn annotate_atoms(&mut self, index: &ComponentIndex, atoms: &[AtomRef<'_>]) {
        let component = index.component();
        for atom in atoms {
            let Some(name) = atom.name() else {
                continue;
            };
            let Some(expected) = index.atom(name) else {
                continue;
            };
            let Some(slot) = self.chemistry.get_mut(atom.index().as_usize()) else {
                continue;
            };
            *slot = Some(AtomChemistry {
                kind: component.kind,
                aromatic: expected.aromatic,
                charge: expected.charge,
                donor: index.is_donor(name),
                acceptor: index.is_acceptor(name),
                weak_acceptor: index.is_weak_acceptor(name),
                stereo: observed_stereo(component, atoms, *atom, expected.stereo),
            });
        }
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
    // Scanning the residue's atoms once per bond endpoint costs
    // `O(bonds x atoms^2)` name comparisons. Ordering the observed atoms by
    // name once turns each endpoint into a binary search, and the sort is
    // stable so atoms sharing a name keep their original relative order.
    let mut by_name: Vec<(&str, AtomRef<'_>)> = atoms
        .iter()
        .filter_map(|atom| atom.name().map(|name| (name, *atom)))
        .collect();
    by_name.sort_by_key(|(name, _)| *name);

    for bond in component.bonds.iter() {
        for (_, atom_a) in named(&by_name, &bond.atom_a) {
            for (_, atom_b) in named(&by_name, &bond.atom_b) {
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

/// The contiguous run of name-ordered atoms called `name`.
fn named<'a, 'b>(
    by_name: &'a [(&'a str, AtomRef<'b>)],
    name: &str,
) -> &'a [(&'a str, AtomRef<'b>)] {
    let start = by_name.partition_point(|(observed, _)| *observed < name);
    let end = by_name.partition_point(|(observed, _)| *observed <= name);
    match by_name.get(start..end) {
        Some(run) => run,
        None => &[],
    }
}

fn alt_compatible(atom_a: AtomRef<'_>, atom_b: AtomRef<'_>) -> bool {
    match (atom_a.alt_id(), atom_b.alt_id()) {
        (Some(a), Some(b)) => a.is_blank() || b.is_blank() || a == b,
        _ => true,
    }
}

/// Classifies a chain by the polymer kind most of its polymer residues share.
///
/// Waters, ions, ligands and other non-polymer components never vote, so a
/// protein chain carrying a few modified or hetero residues is still a protein.
/// A chain is left unclassified when no polymer kind holds a strict majority
/// of the polymer residues.
fn classify_chain(
    data: &mut molframe_core::StructureData,
    chain: ChainIndex,
    kinds: &BTreeMap<ComponentKind, usize>,
) -> Result<(), Diagnostic> {
    let polymer = |kind: ComponentKind| match kinds.get(&kind) {
        Some(count) => *count,
        None => 0,
    };
    let amino = polymer(ComponentKind::AminoAcid);
    let nucleic = polymer(ComponentKind::Nucleotide);
    let sugar = polymer(ComponentKind::Saccharide);
    let total = amino + nucleic + sugar;
    let kind = if amino * 2 > total {
        PolymerKind::Protein
    } else if nucleic * 2 > total {
        PolymerKind::NucleicHybrid
    } else if sugar * 2 > total {
        PolymerKind::Saccharide
    } else {
        return Ok(());
    };
    data.topology
        .chains
        .set_polymer_kind(chain, kind)
        .map_err(|error| Diagnostic::new(Code::E3001).with_context("cause", error.to_string()))
}

#[cfg(test)]
#[path = "annotate_tests.rs"]
mod tests;

//! Applying component chemistry to an immutable structure snapshot.

mod atom_chemistry;
mod component_bonds;
mod component_index;
mod hbond_roles;
mod linkage;
mod resolver;
mod tally;
mod worker;

use component_bonds::alt_compatible;
use linkage::link_polymer;
use rayon::prelude::*;
use resolver::Resolver;
use tally::Tally;
use worker::{UnitOutput, Worker};

use crate::{ComponentKind, ComponentProvider};
pub use atom_chemistry::HBOND_WEAK_ACCEPTOR_ANNOTATION;
use atom_chemistry::{AtomChemistry, attach_annotations};
use molframe_core::bond::BondTableBuilder;
use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::{Code, Diagnostic, Diagnostics};
use molframe_core::index::ChainIndex;
use molframe_core::structure::{ResidueRef, Structure};
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

/// Residues handled by one parallel work unit.
///
/// Units are large enough to amortise scheduling and small enough that one huge
/// chain still spreads over every core.
const UNIT_RESIDUES: usize = 256;

/// Adds CCD internal bonds and checked polymer links without replacing file
/// connectivity.
///
/// Residues are independent, so they are processed in parallel units. Units are
/// merged in file order, which keeps the bond table, annotations and findings
/// identical to a sequential pass. Findings that repeat per residue (`W3202`,
/// `W3301`) are reported once each, with exact counts and a bounded sample.
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
    let chains: Vec<(ChainIndex, Vec<ResidueRef<'_>>)> = structure
        .data()
        .chains()
        .map(|chain| (chain.index(), chain.residues().collect()))
        .collect();
    let mut units: Vec<&[ResidueRef<'_>]> = Vec::new();
    let mut chain_units = Vec::with_capacity(chains.len());
    for (_, residues) in &chains {
        let first = units.len();
        units.extend(residues.chunks(UNIT_RESIDUES));
        chain_units.push(first..units.len());
    }
    let resolver = Resolver::new(provider);
    let outputs: Vec<Result<UnitOutput, Diagnostic>> = units
        .par_iter()
        .map_init(Worker::default, |worker, residues| {
            worker.run(&resolver, residues)
        })
        .collect();

    let mut outputs = outputs.into_iter();
    let mut chemistry: Vec<Option<AtomChemistry>> = vec![None; structure.atom_count() as usize];
    let mut inconsistent = Tally::default();
    let mut missing_atoms = Tally::default();
    let mut missing = BTreeSet::new();
    let mut findings = Diagnostics::new();
    let components = match &polymer_link_policy {
        PolymerLinkPolicy::Disabled => BTreeMap::new(),
        PolymerLinkPolicy::Explicit { .. } => resolver.components(),
    };
    for ((chain, residues), range) in chains.iter().zip(&chain_units) {
        let mut kinds = BTreeMap::new();
        for output in outputs.by_ref().take(range.len()) {
            let output = output?;
            for (kind, count) in output.kinds {
                *kinds.entry(kind).or_insert(0) += count;
            }
            for (atom, value) in output.chemistry {
                if let Some(slot) = chemistry.get_mut(atom) {
                    *slot = Some(value);
                }
            }
            for bond in output.bonds {
                bonds.push(bond);
            }
            inconsistent.merge(output.inconsistent);
            missing_atoms.merge(output.missing_atoms);
            missing.extend(output.missing);
        }
        classify_chain(&mut data, *chain, &kinds)?;
        link_polymer(
            residues,
            &components,
            &mut bonds,
            &mut findings,
            &polymer_link_policy,
        );
    }
    if let Some(finding) = inconsistent.into_diagnostic(Code::W3202, None) {
        findings.push(finding);
    }
    if let Some(finding) = missing_atoms.into_diagnostic(Code::W3301, Some("missing atoms")) {
        findings.push(finding);
    }
    for component in &missing {
        findings
            .push(Diagnostic::new(Code::W3201).with_context("component", component.to_string()));
    }
    data.bonds = bonds.finish();
    hbond_roles::refine_linked_acceptors(&data, &mut chemistry);
    attach_annotations(&mut data, &chemistry)?;
    for finding in molframe_core::structure::validate(&data) {
        findings.push(finding);
    }
    Ok(ChemistryReport {
        structure: Structure::new(data),
        findings: findings.finish(),
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

//! Structure projection for chemistry-derived side-chain torsions.

use crate::{AnalysisPolicy, Diagnostic, Findings, ResidueIndex, Structure};
use molframe_chem::{ComponentProvider, PolymerAtomRole, SideChainRoles, side_chain_definition};
use molframe_core::Code;
use molframe_core::contract::DictionaryVersion;
use molframe_core::structure::{ResidueRef, Structure as CoreStructure};

use super::roles::{atom_role, require_polymer_roles};

/// χ torsions and their defining component atom path for one residue.
#[derive(Clone, Debug, PartialEq)]
pub struct SideChainTorsionRecord {
    /// Residue position in structure topology.
    pub residue: ResidueIndex,
    /// Ordered CCD-local atom identifiers beginning with the annotated anchors.
    pub atoms: Box<[Box<str>]>,
    /// Consecutive χ1–χ5 values in radians; missing geometry remains `None`.
    pub torsions: Box<[Option<f64>]>,
}

/// Side-chain torsions plus provider/policy traceability.
#[derive(Clone, Debug)]
pub struct SideChainTorsionReport {
    /// Residues for which CCD connectivity defines a side-chain path.
    pub records: Vec<SideChainTorsionRecord>,
    /// Altloc and unknown-component findings.
    pub findings: Vec<Diagnostic>,
    /// Exact component dictionary version used.
    pub dictionary_version: DictionaryVersion,
}

struct ResolvedSideChainRoles<'a> {
    nitrogen: &'a str,
    alpha_carbon: &'a str,
    side_chain_atoms: Vec<&'a str>,
}

/// Computes χ1–χ5 from CCD connectivity and selected structure coordinates.
///
/// # Errors
///
/// Returns a diagnostic when the role annotation is absent or ambiguous, or
/// when the component provider fails. Unknown components are non-fatal findings.
pub fn structure_side_chain_torsions(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    policy: &AnalysisPolicy,
) -> Result<SideChainTorsionReport, Findings> {
    let engine = structure.engine();
    require_polymer_roles(engine)?;
    let selected = engine.resolve_altlocs(policy);
    let mut findings = selected.warnings;
    let mut records = Vec::new();
    for residue in engine.data().residues() {
        let Some(component_id) = residue.name() else {
            continue;
        };
        let Some(component) = provider.get(component_id)? else {
            findings
                .push(Diagnostic::new(crate::Code::W3201).with_context("component", component_id));
            continue;
        };
        let Some(resolved) = side_chain_roles(engine, residue, &selected.value)? else {
            continue;
        };
        let roles = SideChainRoles {
            nitrogen: resolved.nitrogen,
            alpha_carbon: resolved.alpha_carbon,
            side_chain_atoms: &resolved.side_chain_atoms,
        };
        let Some(definition) = side_chain_definition(&component, &roles) else {
            continue;
        };
        let positions = definition
            .atoms
            .iter()
            .map(|name| {
                let mut matches = residue
                    .atoms()
                    .filter(|atom| selected.value.contains(atom.index().get()))
                    .filter(|atom| atom.name() == Some(name));
                let first = matches
                    .next()
                    .and_then(molframe_core::structure::AtomRef::position);
                if matches.next().is_some() {
                    None
                } else {
                    first
                }
            })
            .collect::<Vec<_>>();
        records.push(SideChainTorsionRecord {
            residue: residue.index(),
            atoms: definition.atoms,
            torsions: molframe_geom::path_torsions(&positions, 5).into_boxed_slice(),
        });
    }
    Ok(SideChainTorsionReport {
        records,
        findings,
        dictionary_version: provider.version().clone(),
    })
}

fn side_chain_roles<'a>(
    structure: &'a CoreStructure,
    residue: ResidueRef<'a>,
    selected: &molframe_core::AtomSelection,
) -> Result<Option<ResolvedSideChainRoles<'a>>, Diagnostic> {
    let nitrogen = unique_role_name(
        structure,
        residue,
        selected,
        PolymerAtomRole::PROTEIN_NITROGEN,
    )?;
    let alpha_carbon = unique_role_name(
        structure,
        residue,
        selected,
        PolymerAtomRole::PROTEIN_ALPHA_CARBON,
    )?;
    let side_chain_atoms = residue
        .atoms()
        .filter(|atom| selected.contains(atom.index().get()))
        .filter(|atom| {
            atom_role(structure, *atom)
                .is_some_and(|role| role.intersects(PolymerAtomRole::PROTEIN_SIDECHAIN))
        })
        .filter_map(molframe_core::structure::AtomRef::name)
        .collect::<Vec<_>>();
    match (nitrogen, alpha_carbon, side_chain_atoms.is_empty()) {
        (Some(nitrogen), Some(alpha_carbon), false) => Ok(Some(ResolvedSideChainRoles {
            nitrogen,
            alpha_carbon,
            side_chain_atoms,
        })),
        _ => Ok(None),
    }
}

fn unique_role_name<'a>(
    structure: &'a CoreStructure,
    residue: ResidueRef<'a>,
    selected: &molframe_core::AtomSelection,
    required: PolymerAtomRole,
) -> Result<Option<&'a str>, Diagnostic> {
    let mut matching = residue
        .atoms()
        .filter(|atom| selected.contains(atom.index().get()))
        .filter(|atom| atom_role(structure, *atom).is_some_and(|role| role.intersects(required)));
    let first = matching
        .next()
        .and_then(molframe_core::structure::AtomRef::name);
    if matching.next().is_some() {
        Err(Diagnostic::new(Code::E4002)
            .with_context("residue", residue.index().to_string())
            .with_context(
                "required",
                format!("unique polymer role {}", required.code()),
            ))
    } else {
        Ok(first)
    }
}

#[cfg(test)]
#[path = "side_chains_tests.rs"]
mod tests;

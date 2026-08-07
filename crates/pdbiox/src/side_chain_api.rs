//! Structure projection for chemistry-derived side-chain torsions.

use crate::{AnalysisPolicy, Diagnostic, ResidueIndex, Structure};
use pdbiox_chem::{ComponentProvider, side_chain_definition};
use pdbiox_core::contract::DictionaryVersion;

/// χ torsions and their defining component atom path for one residue.
#[derive(Clone, Debug, PartialEq)]
pub struct SideChainTorsionRecord {
    /// Residue position in structure topology.
    pub residue: ResidueIndex,
    /// Ordered atom names beginning `N`, `CA`, `CB`.
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

/// Computes χ1–χ5 from CCD connectivity and selected structure coordinates.
///
/// # Errors
///
/// Returns a provider diagnostic. Unknown components are non-fatal findings.
pub fn structure_side_chain_torsions(
    structure: &Structure,
    provider: &dyn ComponentProvider,
    policy: &AnalysisPolicy,
) -> Result<SideChainTorsionReport, Diagnostic> {
    let selected = structure.resolve_altlocs(policy);
    let mut findings = selected.warnings;
    let mut records = Vec::new();
    for residue in structure.data().residues() {
        let Some(component_id) = residue.name() else {
            continue;
        };
        let Some(component) = provider.get(component_id)? else {
            findings
                .push(Diagnostic::new(crate::Code::W3201).with_context("component", component_id));
            continue;
        };
        let Some(definition) = side_chain_definition(&component) else {
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
                    .and_then(pdbiox_core::structure::AtomRef::position);
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
            torsions: pdbiox_geom::path_torsions(&positions, 5).into_boxed_slice(),
        });
    }
    Ok(SideChainTorsionReport {
        records,
        findings,
        dictionary_version: provider.version().clone(),
    })
}

#[cfg(test)]
#[path = "side_chain_api_tests.rs"]
mod tests;

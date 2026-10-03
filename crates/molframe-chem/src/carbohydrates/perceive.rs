//! Structure-level saccharide inventory using existing component chemistry.

use super::links::collect_links;
use super::rings::residue_rings;
use super::snfg::snfg_symbol;
use super::types::{CarbohydrateOptions, CarbohydrateReport};
use crate::{ComponentKind, ComponentProvider, PolymerLinkPolicy, apply_component_chemistry};
use molframe_core::{BondOrder, BondTableBuilder, Code, Diagnostic, Structure};
use std::collections::BTreeMap;

/// Perceives complete five/six-membered saccharide rings and their attachments.
///
/// A provider adds internal CCD bonds through `apply_component_chemistry`, with
/// file bonds taking precedence. Without a provider, only input topology is
/// used. Recognized CCD aliases and provider-classified saccharides are searched;
/// unknown sugar stereochemistry remains an unknown symbol, never guessed.
/// Alternate conformers are reported separately and never mixed in one ring.
/// Spatial inference is restricted to vacant donor/acceptor sites within 2 Å,
/// and ambiguous competing sites are left unlinked. Geometry uses the first
/// model. The input snapshot is unchanged.
///
/// # Errors
///
/// Returns the existing diagnostic on provider storage failure.
///
/// ```
/// use molframe_chem::{carbohydrates, CarbohydrateOptions};
/// let structure = molframe_core::Structure::new(molframe_core::StructureData::empty());
/// let report = carbohydrates(&structure, None, CarbohydrateOptions::default())?;
/// assert!(report.monosaccharides.is_empty());
/// # Ok::<(), molframe_core::Diagnostic>(())
/// ```
pub fn carbohydrates(
    structure: &Structure,
    provider: Option<&dyn ComponentProvider>,
    options: CarbohydrateOptions,
) -> Result<CarbohydrateReport, Diagnostic> {
    let enriched;
    let structure = match provider {
        Some(provider) => {
            enriched = apply_component_chemistry(structure, provider, PolymerLinkPolicy::Disabled)?
                .structure;
            &enriched
        }
        None => structure,
    };
    let mut report = CarbohydrateReport {
        monosaccharides: Vec::new(),
        links: Vec::new(),
        terminal_links: Vec::new(),
        incomplete_residues: Vec::new(),
        dictionary_version: provider.map(|provider| provider.version().clone()),
    };
    // Double/aromatic edges do not describe a saturated hemiacetal sugar ring.
    // Keep topology traversal separate from linkage provenance, which uses the
    // original/enriched table unchanged.
    let mut ring_bonds = BondTableBuilder::new();
    for bond in structure.data().bonds.iter() {
        if matches!(bond.order, BondOrder::Single | BondOrder::Unknown) {
            ring_bonds.push(bond);
        }
    }
    let ring_bonds = ring_bonds.finish();
    let graph = ring_bonds.adjacency(structure.atom_count());
    let mut components = BTreeMap::new();
    for residue in structure.data().residues() {
        let Some(name) = residue.name() else {
            continue;
        };
        if !components.contains_key(name) {
            let component = match provider {
                Some(provider) => provider.get(name)?,
                None => None,
            };
            components.insert(name, component);
        }
        let component = components
            .get(name)
            .and_then(|component| component.as_deref());
        let symbol = snfg_symbol(name);
        if symbol.is_none()
            && component.is_none_or(|component| component.kind != ComponentKind::Saccharide)
        {
            continue;
        }
        // A corrupted dense graph is not a chemical ring and must not consume
        // an unbounded amount of traversal work.
        if residue
            .atoms()
            .any(|atom| graph.neighbours(atom.index()).len() > 64)
        {
            return Err(Diagnostic::new(Code::E6008)
                .with_message("saccharide graph exceeds traversal degree limit 64"));
        }
        let start = report.monosaccharides.len();
        residue_rings(
            structure,
            graph,
            residue,
            component,
            symbol,
            &mut report.monosaccharides,
        );
        if report.monosaccharides.len() == start {
            report.incomplete_residues.push(residue.index());
        }
    }
    report
        .monosaccharides
        .sort_by(|a, b| (a.residue, &a.ring_atoms).cmp(&(b.residue, &b.ring_atoms)));
    collect_links(structure, &mut report, options.spatial_fallback);
    Ok(report)
}

#[cfg(test)]
#[path = "perceive_tests.rs"]
mod tests;

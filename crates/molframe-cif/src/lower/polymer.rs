//! The polymer kind a file declares for each of its polymer entities.
//!
//! `_entity_poly.type` is the depositor's statement of what a polymer is. A
//! chain of a polymer entity is only known to be *some* polymer until this is
//! read; without it, a protein cannot be told from a nucleic acid except by
//! inspecting every residue against a chemical dictionary. Reading the
//! declaration is not a guess, so it happens here, at lowering, for every read.

use crate::document::DataBlock;
use molframe_core::diagnostic::{Code, Diagnostic, Diagnostics};
use molframe_core::structure::StructureData;
use molframe_core::topology::PolymerKind;

/// Sets each polymer chain's kind from its entity's declared `_entity_poly.type`.
///
/// Only a chain still marked as an unspecified polymer is changed: a chain the
/// file does not declare a polymer is left alone, and a type this reader does
/// not model stays [`PolymerKind::Other`].
pub(super) fn classify(block: &DataBlock, data: &mut StructureData, findings: &mut Diagnostics) {
    let Some(polymers) = block.category("entity_poly") else {
        return;
    };
    let mut declared = Vec::with_capacity(polymers.row_count());
    for row in 0..polymers.row_count() {
        let (Some(entity_id), Some(kind)) = (
            polymers.identifier("entity_id", row),
            polymers.text("type", row).and_then(polymer_kind),
        ) else {
            continue;
        };
        let Some(symbol) = data.dictionary.get(&entity_id) else {
            continue;
        };
        if let Some(entity) = data.topology.entities.find_by_id(symbol) {
            declared.push((entity, kind));
        }
    }
    if declared.is_empty() {
        return;
    }
    let chains: Vec<_> = data.topology.chains.iter().collect();
    for chain in chains {
        if data.topology.chains.polymer_kind(chain) != Some(PolymerKind::Other) {
            continue;
        }
        let Some(entity) = data.topology.chains.entity(chain) else {
            continue;
        };
        let Some(&(_, kind)) = declared.iter().find(|(declared, _)| *declared == entity) else {
            continue;
        };
        if data.topology.chains.set_polymer_kind(chain, kind).is_err() {
            findings.push(Diagnostic::new(Code::E3001));
        }
    }
}

/// The polymer kind an `_entity_poly.type` value names, where it names one.
///
/// The values are the PDBx/mmCIF dictionary's enumeration. Anything else, or
/// a type with no dedicated kind (`other`, `peptide nucleic acid`,
/// `cyclic-pseudo-peptide`), returns `None` and the chain stays `Other`.
fn polymer_kind(declared: &str) -> Option<PolymerKind> {
    let declared = declared.trim();
    let is = |name: &str| declared.eq_ignore_ascii_case(name);
    if is("polypeptide(L)") || is("polypeptide(D)") {
        Some(PolymerKind::Protein)
    } else if is("polydeoxyribonucleotide") {
        Some(PolymerKind::Dna)
    } else if is("polyribonucleotide") {
        Some(PolymerKind::Rna)
    } else if is("polydeoxyribonucleotide/polyribonucleotide hybrid") {
        Some(PolymerKind::NucleicHybrid)
    } else if is("polysaccharide(D)") || is("polysaccharide(L)") {
        Some(PolymerKind::Saccharide)
    } else {
        None
    }
}

#[cfg(test)]
#[path = "polymer_tests.rs"]
mod tests;

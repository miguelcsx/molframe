//! Component roles for structures read without a component dictionary.
//!
//! A component dictionary is the authority on what a residue is, but most
//! structure files ship without one and every selector that asks "is this
//! protein?" would then answer no. The standard amino acids and nucleotides
//! are recognisable from their names alone, so this pass marks those and
//! leaves every other residue unknown. Solvent, ions and ligands are
//! deliberately not guessed from names or atom counts: the file's declared
//! entities say what they are, and a wrong guess here would silently override
//! them. It never replaces a column a dictionary already wrote.
//!
//! One pass over the residues, `O(N)` in atoms, allocating a single column.

use crate::ComponentKind;
use crate::standard_bonds::{is_amino_acid_component, is_nucleotide_component};
use molframe_core::annotation::{AnnotationColumn, AtomAnnotation};
use molframe_core::column::Presence;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::structure::Structure;

/// The role the built-in tables assign to a residue named `component`:
/// amino acid or nucleotide for the standard components, unknown otherwise.
#[must_use]
pub fn standard_component_kind(component: &str) -> ComponentKind {
    if is_amino_acid_component(component) {
        ComponentKind::AminoAcid
    } else if is_nucleotide_component(component) {
        ComponentKind::Nucleotide
    } else {
        ComponentKind::Unknown
    }
}

/// Fills the component-role column from the built-in tables when absent.
///
/// A structure that already carries the column is returned unchanged, sharing
/// all storage.
///
/// # Errors
///
/// Returns a diagnostic when the column exceeds the supported row count.
pub fn annotate_standard_components(structure: &Structure) -> Result<Structure, Diagnostic> {
    if structure
        .annotations()
        .get(molframe_core::COMPONENT_KIND_ANNOTATION)
        .is_some()
    {
        return Ok(structure.clone());
    }
    let mut entries =
        vec![(ComponentKind::Unknown.code(), Presence::Unknown); structure.atom_count() as usize];
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            let kind = residue
                .name()
                .map_or(ComponentKind::Unknown, standard_component_kind);
            if kind == ComponentKind::Unknown {
                continue;
            }
            for atom in residue.atoms() {
                if let Some(entry) = entries.get_mut(atom.index().as_usize()) {
                    *entry = (kind.code(), Presence::Present);
                }
            }
        }
    }
    let column = AnnotationColumn::from_entries(entries).map_err(|error| {
        Diagnostic::new(Code::E6008)
            .with_message("component-role column exceeds the supported row count")
            .with_context("cause", error.to_string())
    })?;
    let mut data = structure.data().clone();
    let _ = data.annotations.insert(
        molframe_core::COMPONENT_KIND_ANNOTATION,
        AtomAnnotation::Integer(column),
    );
    Ok(Structure::new(data))
}

#[cfg(test)]
#[path = "standard_components_tests.rs"]
mod tests;

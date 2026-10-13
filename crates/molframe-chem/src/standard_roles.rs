//! Built-in protein backbone roles for structures read without a dictionary.
//!
//! A researcher who asks for secondary structure expects the wwPDB atom names
//! to mean what they have always meant: `N`, `CA`, `C` and `O` are the peptide
//! backbone and `CB` is the first side-chain carbon. This pass writes exactly
//! that, and only for residues that are amino acids, either because the
//! structure already carries a component-kind column (a dictionary decided) or
//! because the residue is one of the standard or common amino acids. Everything
//! else keeps an unknown role rather than a guessed one, so a ligand that
//! happens to own an atom called `C` is never read as a carbonyl carbon. A
//! caller who needs other names supplies an explicit
//! [`PolymerRoleProfile`](crate::PolymerRoleProfile) instead.
//!
//! One pass over the residues, `O(N)` in atoms, allocating a single column.

use crate::standard_bonds::is_amino_acid_component;
use crate::{ComponentKind, PolymerAtomRole};
use molframe_core::annotation::{AnnotationColumn, AtomAnnotation};
use molframe_core::column::Presence;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::structure::Structure;

/// Identifier of the built-in protein role assignment.
pub const STANDARD_PROTEIN_ROLE_PROFILE: &str = "standard-protein-atom-names-v1";

/// The role the standard protein profile gives to an atom of an amino acid.
#[must_use]
pub fn standard_protein_atom_role(atom_name: &str) -> PolymerAtomRole {
    match atom_name {
        "N" => PolymerAtomRole::PROTEIN_NITROGEN,
        "CA" => PolymerAtomRole::PROTEIN_ALPHA_CARBON,
        "C" => PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
        "O" => PolymerAtomRole::PROTEIN_CARBONYL_OXYGEN,
        "CB" => PolymerAtomRole::PROTEIN_BETA_CARBON,
        _ => PolymerAtomRole::UNKNOWN,
    }
}

/// Fills the polymer-role column from the standard protein profile when absent.
///
/// A structure that already carries the column is returned unchanged, sharing
/// all storage, so an explicit profile always wins over this default.
///
/// # Errors
///
/// Returns a diagnostic when the column exceeds the supported row count.
pub fn annotate_standard_protein_roles(structure: &Structure) -> Result<Structure, Diagnostic> {
    if structure
        .annotations()
        .get(molframe_core::POLYMER_ATOM_ROLE_ANNOTATION)
        .is_some()
    {
        return Ok(structure.clone());
    }
    let kinds = match structure
        .annotations()
        .get(molframe_core::COMPONENT_KIND_ANNOTATION)
    {
        Some(AtomAnnotation::Integer(column)) => Some(column),
        _ => None,
    };
    let mut entries =
        vec![(PolymerAtomRole::UNKNOWN.code(), Presence::Unknown); structure.atom_count() as usize];
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            let declared = kinds.and_then(|column| {
                let first = residue.atoms().next()?;
                match column.get(first.index().get()) {
                    Some((code, Presence::Present)) => ComponentKind::from_code(code),
                    _ => None,
                }
            });
            let amino = match declared {
                Some(kind) => kind == ComponentKind::AminoAcid,
                None => residue.name().is_some_and(is_amino_acid_component),
            };
            if !amino {
                continue;
            }
            for atom in residue.atoms() {
                let role = match atom.name() {
                    Some(name) => standard_protein_atom_role(name),
                    None => PolymerAtomRole::UNKNOWN,
                };
                let presence = if role == PolymerAtomRole::UNKNOWN {
                    Presence::Inapplicable
                } else {
                    Presence::Present
                };
                if let Some(entry) = entries.get_mut(atom.index().as_usize()) {
                    *entry = (role.code(), presence);
                }
            }
        }
    }
    let column = AnnotationColumn::from_entries(entries).map_err(|error| {
        Diagnostic::new(Code::E6009)
            .with_message("polymer-role column exceeds the supported row count")
            .with_context("cause", error.to_string())
    })?;
    let mut data = structure.data().clone();
    let _ = data.annotations.insert(
        molframe_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(column),
    );
    Ok(Structure::new(data))
}

#[cfg(test)]
#[path = "standard_roles_tests.rs"]
mod tests;

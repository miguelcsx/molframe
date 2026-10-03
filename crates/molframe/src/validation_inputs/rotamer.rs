//! Side-chain torsion definitions that name their reference distributions.

use super::ValidationInputError;
use crate::document::read_document;
use molframe_validate::{RotamerDefinition, RotamerProfile};
use serde::Deserialize;
use std::path::Path;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Document {
    profile: Identity,
    definition: Vec<Definition>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Identity {
    id: Box<str>,
    version: Box<str>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Definition {
    component_id: Box<str>,
    chi_index: u8,
    atoms: [Box<str>; 4],
    distribution: Box<str>,
}

/// Reads a rotamer profile from JSON or TOML.
///
/// A `[profile]` table names the profile (`id`, `version`), and each
/// `[[definition]]` gives one component-local torsion: the CCD `component_id`,
/// the one-based `chi_index`, the four component `atoms` that define the
/// dihedral, and the `distribution` it is assessed against in the reference
/// library the caller pairs it with. Keys equal the fields of
/// [`RotamerDefinition`]; unknown keys are errors.
///
/// ```toml
/// [profile]
/// id = "chi-wells"
/// version = "2026-10"
///
/// [[definition]]
/// component_id = "SER"
/// chi_index = 1
/// atoms = ["N", "CA", "CB", "OG"]
/// distribution = "ser-chi1"
/// ```
///
/// # Errors
///
/// Returns [`ValidationInputError`] for I/O, syntax or schema errors, or when
/// [`RotamerProfile::new`] refuses the profile (empty identity, a zero torsion
/// index, a duplicate component and torsion pair).
pub fn read_rotamer_profile(
    path: impl AsRef<Path>,
) -> Result<RotamerProfile, ValidationInputError> {
    let document: Document = read_document(path.as_ref())?;
    let definitions = document
        .definition
        .into_iter()
        .map(|definition| RotamerDefinition {
            component_id: definition.component_id,
            chi_index: definition.chi_index,
            atoms: definition.atoms,
            distribution: definition.distribution,
        });
    Ok(RotamerProfile::new(
        document.profile.id,
        document.profile.version,
        definitions,
    )?)
}

#[cfg(test)]
#[path = "rotamer_tests.rs"]
mod tests;

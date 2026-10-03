//! Restraints and models that govern atoms chosen by selection expressions.

use super::ValidationInputError;
use crate::AnalysisPolicy;
use crate::document::read_document;
use crate::structure::{QueryStructure, Structure};
use molframe_core::AtomSelection;

use molframe_validate::{PlaneRestraint, TlsGroup, TlsModel};
use serde::Deserialize;
use std::collections::BTreeSet;
use std::path::Path;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Planes {
    plane: Vec<Plane>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Plane {
    id: String,
    selection: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Groups {
    group: Vec<Group>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Group {
    id: String,
    selection: String,
    origin: [f64; 3],
    t: [[f64; 3]; 3],
    l: [[f64; 3]; 3],
    s: [[f64; 3]; 3],
}

/// Evaluates one selection expression, refusing an empty result.
fn select(
    structure: &Structure,
    policy: &AnalysisPolicy,
    id: &str,
    source: &str,
) -> Result<AtomSelection, ValidationInputError> {
    let evaluation = structure
        .engine()
        .select_text(source, policy)
        .map_err(|findings| ValidationInputError::Selection {
            id: id.to_owned(),
            findings,
        })?;
    if evaluation.selection.is_empty() {
        return Err(ValidationInputError::invalid(
            format!("selection for {id}"),
            "the expression selects no atoms in this structure",
        ));
    }
    Ok(evaluation.selection)
}

/// Refuses an empty or repeated identifier, so a report can name its subject.
fn claim_id<'a>(
    seen: &mut BTreeSet<&'a str>,
    kind: &str,
    id: &'a str,
) -> Result<(), ValidationInputError> {
    if id.trim().is_empty() {
        return Err(ValidationInputError::invalid(
            format!("{kind}.id"),
            "must not be empty",
        ));
    }
    if !seen.insert(id) {
        return Err(ValidationInputError::invalid(
            format!("{kind}.id"),
            format!("`{id}` is used more than once"),
        ));
    }
    Ok(())
}

/// Reads explicit plane restraints and resolves their atoms in `structure`.
///
/// Each `[[plane]]` has a unique `id` and a `selection` expression in the query
/// language, evaluated under `policy` (which decides, for instance, whether
/// `chain A` reads label or author identifiers). A selection that matches no
/// atom is an error, because a restraint that silently covers nothing passes
/// every structure.
///
/// ```toml
/// [[plane]]
/// id = "peptide-1"
/// selection = "resid 1:2 and name CA C O N"
/// ```
///
/// # Errors
///
/// Returns [`ValidationInputError`] for I/O, syntax or schema errors, duplicate
/// ids, or a selection that fails or is empty.
pub fn read_plane_restraints(
    path: impl AsRef<Path>,
    structure: &Structure,
    policy: &AnalysisPolicy,
) -> Result<Vec<PlaneRestraint>, ValidationInputError> {
    let document: Planes = read_document(path.as_ref())?;
    let mut seen = BTreeSet::new();
    let mut restraints = Vec::with_capacity(document.plane.len());
    for plane in &document.plane {
        claim_id(&mut seen, "plane", &plane.id)?;
        restraints.push(PlaneRestraint {
            id: plane.id.clone(),
            atoms: select(structure, policy, &plane.id, &plane.selection)?,
        });
    }
    Ok(restraints)
}

/// Reads TLS groups and resolves their atoms in `structure`.
///
/// Each `[[group]]` has a unique `id`, a `selection` expression evaluated under
/// `policy`, the TLS `origin` (three coordinates in the structure's units) and
/// the three tensors `t`, `l` and `s`, each as three rows of three numbers, in
/// the crystallographic TLS convention the validation kernel documents. Only
/// shape and finiteness are checked here; the kernel decides whether the
/// tensors are admissible.
///
/// ```toml
/// [[group]]
/// id = "A-1"
/// selection = "chain A"
/// origin = [0.0, 0.0, 0.0]
/// t = [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]]
/// l = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
/// s = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
/// ```
///
/// # Errors
///
/// Returns [`ValidationInputError`] for I/O, syntax or schema errors, duplicate
/// ids, a non-finite tensor entry, or a selection that fails or is empty.
pub fn read_tls_groups(
    path: impl AsRef<Path>,
    structure: &Structure,
    policy: &AnalysisPolicy,
) -> Result<Vec<TlsGroup>, ValidationInputError> {
    let document: Groups = read_document(path.as_ref())?;
    let mut seen = BTreeSet::new();
    let mut groups = Vec::with_capacity(document.group.len());
    for group in &document.group {
        claim_id(&mut seen, "group", &group.id)?;
        let finite = group
            .origin
            .iter()
            .chain(group.t.iter().chain(&group.l).chain(&group.s).flatten())
            .all(|value| value.is_finite());
        if !finite {
            return Err(ValidationInputError::invalid(
                format!("group {}", group.id),
                "origin and tensors must be finite",
            ));
        }
        groups.push(TlsGroup {
            id: group.id.clone(),
            atoms: select(structure, policy, &group.id, &group.selection)?,
            model: TlsModel {
                origin: group.origin,
                translation: group.t,
                libration: group.l,
                screw: group.s,
            },
        });
    }
    Ok(groups)
}

#[cfg(test)]
#[path = "selections_tests.rs"]
mod tests;

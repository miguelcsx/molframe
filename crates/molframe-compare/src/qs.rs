//! The published QS-score for the interface between two chains.
//!
//! This follows the quaternary-structure score of Bertoni et al. (2017), as
//! implemented in `OpenStructure`. Each residue is represented by its CB atom, or
//! its CA atom for glycine; a residue without either takes no part. Two residues
//! of different chains are in contact when their representative atoms lie closer
//! than the contact distance (12 Å). A contact at distance `d` has weight
//!
//! ```text
//! w(d) = 1                          for d <= 5
//! w(d) = exp(-2 ((d - 5) / 4.28)^2)  for d > 5
//! ```
//!
//! and the score compares the reference (native) contact set `R`, with distances
//! `d_ref`, with the model set `M`, with distances `d_mod`:
//!
//! ```text
//!                 sum over shared (i,j) of  w(d_ref) * w(d_mod)
//! QS = ------------------------------------------------------------------
//!      (sum over R of w(d_ref)^2  +  sum over M of w(d_mod)^2) / 2
//! ```
//!
//! A contact is shared when the same residue pair is in contact in both. An
//! identical interface scores 1 (numerator and denominator are both the sum of
//! squared weights), one with no shared contact scores 0, and moving a shared
//! contact away lowers its weight product and so the score, even when it is the
//! only contact. Contacts present on one side only enter the denominator alone.
//!
//! The contact definition and the weight function are those of the published
//! score. The normalisation is the point to treat with care: it is chosen here as
//! the mean of the two squared-weight masses, which has the properties above, and
//! it has not been checked numerically against `OpenStructure` output.
//!
//! Model and native share atom numbering, so a residue is the same in both. The
//! atom-contact Jaccard overlap this crate used to call the QS-score lives in
//! [`atom_contact_jaccard`] and is not comparable with published values.

use std::collections::BTreeMap;

use molframe_core::contract::Namespace;
use molframe_core::structure::Structure;

use crate::CompareError;
use crate::interface::{AtomInfo, UNANNOTATED, atom_infos, chain_atoms};

#[path = "qs_atom_contact.rs"]
mod atom_contact;

pub use atom_contact::{
    AtomContactJaccardOptions, atom_contact_jaccard, atom_contact_jaccard_in_namespace,
};

/// Distance up to which a contact has full weight, in ångström.
const WEIGHT_PLATEAU: f64 = 5.0;
/// Width of the Gaussian fall-off of the contact weight, in ångström.
const WEIGHT_WIDTH: f64 = 4.28;

/// Result policy when neither structure has an interface contact.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EmptyQsPolicy {
    /// Two empty contact sets agree perfectly.
    Perfect,
    /// Treat an empty interface domain as an error.
    Error,
}

/// Explicit definition of a quaternary-structure comparison.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QsOptions {
    /// Distance between residue representative atoms (CB, or CA for glycine)
    /// below which two residues of different chains are in contact, in ångström.
    /// The published value is 12.
    pub contact_distance: f32,
    /// Behavior when both contact sets are empty.
    pub empty_policy: EmptyQsPolicy,
}

impl QsOptions {
    /// The published definition: contacts below 12 Å.
    #[must_use]
    pub const fn published() -> Self {
        Self::standard(12.0)
    }

    /// Two empty contact sets agree perfectly, at a caller-selected cutoff.
    #[must_use]
    pub const fn standard(contact_distance: f32) -> Self {
        Self {
            contact_distance,
            empty_policy: EmptyQsPolicy::Perfect,
        }
    }
}

/// Scores the interface between two chains in `model` against `native`.
///
/// See the module documentation for the exact formula. Cost is quadratic in the
/// number of residues of the two chains.
///
/// # Errors
///
/// Returns an error for unequal structures, an invalid contact distance, a chain
/// with no residue representative atoms, or an empty domain when
/// [`EmptyQsPolicy::Error`] is selected.
pub fn qs_score(
    model: &Structure,
    native: &Structure,
    first_chain: &str,
    second_chain: &str,
    options: QsOptions,
) -> Result<f64, CompareError> {
    qs_score_in_namespace(
        model,
        native,
        first_chain,
        second_chain,
        Namespace::Label,
        options,
    )
}

/// Scores an interface with chain names interpreted in one namespace.
///
/// # Errors
///
/// Returns the same errors as [`qs_score`] and rejects an explicit namespace.
pub fn qs_score_in_namespace(
    model: &Structure,
    native: &Structure,
    first_chain: &str,
    second_chain: &str,
    namespace: Namespace,
    options: QsOptions,
) -> Result<f64, CompareError> {
    if model.atom_count() != native.atom_count() {
        return Err(CompareError::LengthMismatch {
            model: model.atom_count() as usize,
            reference: native.atom_count() as usize,
        });
    }
    let first = chain_atoms(native, first_chain, namespace)?;
    let second = chain_atoms(native, second_chain, namespace)?;
    qs_on_positions(
        model.positions(),
        native.positions(),
        &first,
        &second,
        &atom_infos(native),
        options,
    )
}

/// The contact weight at a representative-atom distance.
pub(crate) fn contact_weight(distance: f64) -> f64 {
    if distance <= WEIGHT_PLATEAU {
        return 1.0;
    }
    let scaled = (distance - WEIGHT_PLATEAU) / WEIGHT_WIDTH;
    (-2.0 * scaled * scaled).exp()
}

/// Scores contact overlap between two coordinate sets, row for row.
///
/// `first` and `second` index rows of both sets and `info` annotates each row.
///
/// # Errors
///
/// Returns an error for an invalid contact distance, a chain without any
/// representative atom, or an empty domain when [`EmptyQsPolicy::Error`] is
/// selected.
pub(crate) fn qs_on_positions(
    model: &[[f32; 3]],
    native: &[[f32; 3]],
    first: &[usize],
    second: &[usize],
    info: &[AtomInfo],
    options: QsOptions,
) -> Result<f64, CompareError> {
    if !options.contact_distance.is_finite() || options.contact_distance <= 0.0 {
        return Err(CompareError::InvalidDistanceCutoff);
    }
    let first = representatives(first, info);
    let second = representatives(second, info);
    if first.is_empty() || second.is_empty() {
        // Without representative atoms the interface cannot be measured, which
        // is different from measuring it and finding no contact.
        return Err(CompareError::NoComparablePairs);
    }
    let reference = contact_distances(native, &first, &second, options.contact_distance);
    let modelled = contact_distances(model, &first, &second, options.contact_distance);

    if reference.is_empty() && modelled.is_empty() {
        return match options.empty_policy {
            EmptyQsPolicy::Perfect => Ok(1.0),
            EmptyQsPolicy::Error => Err(CompareError::NoComparablePairs),
        };
    }
    let reference_mass: f64 = reference.values().map(|&d| contact_weight(d).powi(2)).sum();
    let model_mass: f64 = modelled.values().map(|&d| contact_weight(d).powi(2)).sum();
    let shared: f64 = reference
        .iter()
        .filter_map(|(pair, &d_ref)| {
            modelled
                .get(pair)
                .map(|&d_mod| contact_weight(d_ref) * contact_weight(d_mod))
        })
        .sum();
    let denominator = f64::midpoint(reference_mass, model_mass);
    if denominator > 0.0 {
        Ok(shared / denominator)
    } else {
        Ok(0.0)
    }
}

/// `(residue, row)` of each residue's representative atom among `rows`.
fn representatives(rows: &[usize], info: &[AtomInfo]) -> Vec<(usize, usize)> {
    rows.iter()
        .filter_map(|&row| {
            let annotation = match info.get(row) {
                Some(&found) => found,
                None => UNANNOTATED,
            };
            annotation
                .representative
                .then_some((annotation.residue, row))
        })
        .collect()
}

/// Distance of every residue pair, keyed by the pair, below `cutoff`.
fn contact_distances(
    positions: &[[f32; 3]],
    first: &[(usize, usize)],
    second: &[(usize, usize)],
    cutoff: f32,
) -> BTreeMap<(usize, usize), f64> {
    let limit = f64::from(cutoff) * f64::from(cutoff);
    let mut pairs = BTreeMap::new();
    for &(a, row_a) in first {
        let Some(&point_a) = positions.get(row_a) else {
            continue;
        };
        for &(b, row_b) in second {
            let Some(&point_b) = positions.get(row_b) else {
                continue;
            };
            let squared = molframe_geom::distance_squared(point_a, point_b);
            if a != b && squared < limit {
                pairs.insert((a.min(b), a.max(b)), squared.sqrt());
            }
        }
    }
    pairs
}

#[cfg(test)]
#[path = "qs_tests.rs"]
mod tests;

//! The quaternary-structure score for an interface between two chains.
//!
//! QS-score asks how much of a native interface a model reproduces, counted in
//! contacts: the atom pairs across the two chains that are in contact in the
//! native, and how many of them the model keeps, against the union of both sets.
//! It is the Jaccard overlap of the two contact sets, so an identical model
//! scores 1 and a model that shares no interface contacts scores 0.
//!
//! Model and native share atom numbering — the model is a rearrangement of the
//! same atoms — so a contact is the same index pair in both.

use std::collections::BTreeSet;

use molframe_core::contract::Namespace;
use molframe_core::structure::Structure;

use crate::CompareError;
use crate::interface::{chain_atoms, contacts};
use crate::numeric::usize_to_f64;

/// Result policy when neither structure has an interface contact.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EmptyQsPolicy {
    /// Standard Jaccard convention: two empty contact sets agree perfectly.
    Perfect,
    /// Treat an empty interface domain as an error.
    Error,
}

/// Explicit definition of a quaternary-structure contact comparison.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QsOptions {
    /// Maximum inter-chain contact distance in ångström.
    pub contact_distance: f32,
    /// Behavior when both contact sets are empty.
    pub empty_policy: EmptyQsPolicy,
}

impl QsOptions {
    /// Standard QS/Jaccard empty-set convention at a caller-selected cutoff.
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
/// `cutoff` is the contact distance, about 5 Å by convention. When neither
/// structure has any interface contact the score is `1.0`, since there is
/// nothing to disagree about.
///
/// Runs in `O(chain a · chain b)` time.
///
/// # Errors
///
/// Returns an error for unequal structures, an invalid contact distance, or an
/// empty domain when [`EmptyQsPolicy::Error`] is selected.
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
    if !options.contact_distance.is_finite() || options.contact_distance <= 0.0 {
        return Err(CompareError::InvalidDistanceCutoff);
    }
    let first = chain_atoms(native, first_chain, namespace)?;
    let second = chain_atoms(native, second_chain, namespace)?;

    let native_contacts: BTreeSet<(usize, usize)> = contacts(
        native.positions(),
        &first,
        &second,
        options.contact_distance,
    )
    .into_iter()
    .collect();
    let model_contacts: BTreeSet<(usize, usize)> =
        contacts(model.positions(), &first, &second, options.contact_distance)
            .into_iter()
            .collect();

    let shared = native_contacts.intersection(&model_contacts).count();
    let union = native_contacts.union(&model_contacts).count();
    if union == 0 {
        match options.empty_policy {
            EmptyQsPolicy::Perfect => Ok(1.0),
            EmptyQsPolicy::Error => Err(CompareError::NoComparablePairs),
        }
    } else {
        Ok(usize_to_f64(shared) / usize_to_f64(union))
    }
}

#[cfg(test)]
#[path = "qs_tests.rs"]
mod tests;

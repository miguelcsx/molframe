//! The atom-contact Jaccard overlap of an interface between two chains.
//!
//! This is the contact-set overlap that earlier versions called the QS-score. It
//! counts atom pairs across the two chains within a distance and takes the
//! Jaccard overlap of the native and model sets. It is not the published
//! QS-score, which works on residue contacts with distance-dependent weights;
//! use [`crate::qs_score`] for that.
//!
//! Model and native share atom numbering, so a contact is the same index pair in
//! both.

use std::collections::BTreeSet;

use molframe_core::contract::Namespace;
use molframe_core::structure::Structure;

use crate::CompareError;
use crate::interface::{chain_atoms, contacts};
use crate::numeric::usize_to_f64;
use crate::quaternary::EmptyQsPolicy;

/// Explicit definition of an atom-contact Jaccard comparison.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomContactJaccardOptions {
    /// Maximum inter-chain contact distance in ångström.
    pub contact_distance: f32,
    /// Behavior when both contact sets are empty.
    pub empty_policy: EmptyQsPolicy,
}

impl AtomContactJaccardOptions {
    /// Standard Jaccard empty-set convention at a caller-selected cutoff.
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
pub fn atom_contact_jaccard(
    model: &Structure,
    native: &Structure,
    first_chain: &str,
    second_chain: &str,
    options: AtomContactJaccardOptions,
) -> Result<f64, CompareError> {
    atom_contact_jaccard_in_namespace(
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
/// Returns the same errors as [`atom_contact_jaccard`] and rejects an explicit namespace.
pub fn atom_contact_jaccard_in_namespace(
    model: &Structure,
    native: &Structure,
    first_chain: &str,
    second_chain: &str,
    namespace: Namespace,
    options: AtomContactJaccardOptions,
) -> Result<f64, CompareError> {
    if model.atom_count() != native.atom_count() {
        return Err(CompareError::LengthMismatch {
            model: model.atom_count() as usize,
            reference: native.atom_count() as usize,
        });
    }
    let first = chain_atoms(native, first_chain, namespace)?;
    let second = chain_atoms(native, second_chain, namespace)?;
    atom_contact_jaccard_on_positions(
        model.positions(),
        native.positions(),
        &first,
        &second,
        options,
    )
}

/// Scores contact overlap between two coordinate sets, row for row.
///
/// `first` and `second` index rows of both sets.
///
/// # Errors
///
/// Returns an error for an invalid contact distance, or for an empty domain
/// when [`EmptyQsPolicy::Error`] is selected.
pub(crate) fn atom_contact_jaccard_on_positions(
    model: &[[f32; 3]],
    native: &[[f32; 3]],
    first: &[usize],
    second: &[usize],
    options: AtomContactJaccardOptions,
) -> Result<f64, CompareError> {
    if !options.contact_distance.is_finite() || options.contact_distance <= 0.0 {
        return Err(CompareError::InvalidDistanceCutoff);
    }
    let native_contacts: BTreeSet<(usize, usize)> =
        contacts(native, first, second, options.contact_distance)
            .into_iter()
            .collect();
    let model_contacts: BTreeSet<(usize, usize)> =
        contacts(model, first, second, options.contact_distance)
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
#[path = "qs_atom_contact_tests.rs"]
mod tests;

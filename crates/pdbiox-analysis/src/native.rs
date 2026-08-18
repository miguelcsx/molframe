//! The fraction of a reference structure's contacts kept by another.
//!
//! Q measures how much of a native contact map survives in a second structure
//! that shares the reference's atom numbering: every atom pair in contact in the
//! reference is checked in the target, and Q is the share still close enough to
//! count. It is the standard order parameter for how folded, or how native-like,
//! a conformation is.
//!
//! A contact is judged kept when the target distance is within a tolerance
//! factor of the reference cutoff, so a modest expansion still counts while a
//! broken contact does not.

use pdbiox_core::structure::Structure;
use pdbiox_spatial::{SpatialBackend, SpatialError, pairs_within_unsorted};

use crate::numeric::usize_to_f64;

/// Why the native-contact fraction could not be computed.
#[derive(Debug, thiserror::Error)]
pub enum NativeError {
    /// The two structures do not share an atom count, so their indices cannot
    /// name the same atoms.
    #[error("reference has {reference} atoms but target has {target}")]
    AtomCountMismatch {
        /// Atom count of the reference.
        reference: u32,
        /// Atom count of the target.
        target: u32,
    },
    /// The neighbour search over the reference failed.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
}

/// How many native contacts there were, how many survived, and their ratio.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NativeContacts {
    /// Contacts present in the reference.
    pub native: usize,
    /// Reference contacts still within range in the target.
    pub kept: usize,
    /// `kept / native`, or `1.0` when the reference has no contacts.
    pub fraction: f64,
}

/// Computes Q between a reference and a target that share atom numbering.
///
/// `cutoff` defines a reference contact; a reference contact counts as kept when
/// the target distance is at most `tolerance * cutoff`. A `tolerance` of `1.0`
/// demands the contact be at least as close as the cutoff again.
///
/// Runs in the cost of one reference contact search plus a pass over its pairs.
///
/// # Errors
///
/// Returns [`NativeError::AtomCountMismatch`] when the structures disagree on
/// atom count, and [`NativeError::Spatial`] from the neighbour search.
pub fn native_contact_fraction(
    reference: &Structure,
    target: &Structure,
    cutoff: f32,
    tolerance: f32,
    backend: SpatialBackend,
) -> Result<NativeContacts, NativeError> {
    if reference.atom_count() != target.atom_count() {
        return Err(NativeError::AtomCountMismatch {
            reference: reference.atom_count(),
            target: target.atom_count(),
        });
    }

    let all =
        pdbiox_core::selection::AtomSelection::from_sorted((0..reference.atom_count()).collect());
    let native_contacts =
        pairs_within_unsorted(reference.positions(), &all, &all, cutoff, backend, None)?;
    let native = native_contacts.len();
    if native == 0 {
        return Ok(NativeContacts {
            native: 0,
            kept: 0,
            fraction: 1.0,
        });
    }

    let target_positions = target.positions();
    let kept_cutoff_squared = f64::from(tolerance * cutoff).powi(2);
    let mut kept = 0usize;
    for pair in native_contacts {
        let (Some(&a), Some(&b)) = (
            target_positions.get(pair.first as usize),
            target_positions.get(pair.second as usize),
        ) else {
            continue;
        };
        if pdbiox_geom::distance_squared(a, b) <= kept_cutoff_squared {
            kept += 1;
        }
    }

    Ok(NativeContacts {
        native,
        kept,
        fraction: usize_to_f64(kept) / usize_to_f64(native),
    })
}

#[cfg(test)]
#[path = "native_tests.rs"]
mod tests;

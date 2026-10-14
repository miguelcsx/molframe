//! The fraction of a reference structure's contacts kept by another.
//!
//! Q measures how much of a native contact map survives in a second structure:
//! every atom pair in contact in the reference is checked in the target, and Q
//! is the share still close enough to count. It is the standard order parameter
//! for how folded, or how native-like, a conformation is.
//!
//! A contact is native when the reference pair is within `cutoff`, and kept when
//! the same pair in the target is within `cutoff + tolerance`, with `tolerance`
//! an additive distance in ångström. A tolerance of zero therefore keeps exactly
//! the pairs that are still in contact, so a structure compared with itself has
//! `Q = 1`. The contact is a hard cutoff over every atom pair, not the smoothed
//! switching function of the Best-Hummer-Eaton formulation.
//!
//! The target may number its atoms differently from the reference. A
//! correspondence names, for each reference atom, the target atom standing for
//! it; reference atoms without a counterpart take no part in any contact.
use molframe_core::{ExecutionContext, selection::AtomSelection, structure::Structure};
use molframe_spatial::{
    PairQuery, SpatialBackend, SpatialError, SpatialSearchOptions, reduce_pairs_within_unsorted,
};

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
    /// A correspondence does not cover the reference's atoms.
    #[error("correspondence covers {given} atoms but the reference has {reference}")]
    CorrespondenceLength {
        /// Entries in the correspondence.
        given: usize,
        /// Atom count of the reference.
        reference: u32,
    },
    /// A correspondence names a target atom the target does not have.
    #[error("correspondence names target atom {index}, but the target has {target} atoms")]
    TargetAtomAbsent {
        /// The offending target atom index.
        index: u32,
        /// Atom count of the target.
        target: usize,
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
    /// Reference atoms that have a counterpart in the target.
    pub matched_atoms: usize,
    /// `kept / native`, or `1.0` when the reference has no contacts.
    pub fraction: f64,
}

/// Slack on the kept test for the single-precision rounding of the pair search,
/// so a pair the reference search accepted at its cutoff is not lost to the
/// wider arithmetic used here when the target is the reference itself.
const DISTANCE_ROUNDING: f64 = 1e-6;

/// Computes Q between a reference and a target that share atom numbering.
///
/// `cutoff` defines a reference contact; a reference contact counts as kept when
/// the target distance is at most `cutoff + tolerance`, both in ångström.
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
    context: &ExecutionContext,
) -> Result<NativeContacts, NativeError> {
    if reference.atom_count() != target.atom_count() {
        return Err(NativeError::AtomCountMismatch {
            reference: reference.atom_count(),
            target: target.atom_count(),
        });
    }
    let all = AtomSelection::All(reference.atom_count());
    let matched = reference.atom_count() as usize;
    count_kept(
        reference,
        &all,
        matched,
        Some,
        target.positions(),
        (cutoff, tolerance),
        backend,
        context,
    )
}

/// Computes Q when the target numbers its atoms differently.
///
/// `target_of[i]` is the target atom standing for reference atom `i`, or `None`
/// when it has no counterpart; only contacts between matched atoms are counted.
/// `target_positions` are the target coordinates those indices refer to.
///
/// # Errors
///
/// Returns [`NativeError::CorrespondenceLength`] when `target_of` does not hold
/// one entry per reference atom, [`NativeError::TargetAtomAbsent`] when an entry
/// is outside `target_positions`, and [`NativeError::Spatial`] from the search.
pub fn native_contact_fraction_matched(
    reference: &Structure,
    target_positions: &[[f32; 3]],
    target_of: &[Option<u32>],
    cutoff: f32,
    tolerance: f32,
    backend: SpatialBackend,
    context: &ExecutionContext,
) -> Result<NativeContacts, NativeError> {
    if target_of.len() != reference.atom_count() as usize {
        return Err(NativeError::CorrespondenceLength {
            given: target_of.len(),
            reference: reference.atom_count(),
        });
    }
    if let Some(index) = target_of
        .iter()
        .flatten()
        .copied()
        .find(|index| *index as usize >= target_positions.len())
    {
        return Err(NativeError::TargetAtomAbsent {
            index,
            target: target_positions.len(),
        });
    }
    let matched: Vec<u32> = target_of
        .iter()
        .enumerate()
        .filter(|(_, target)| target.is_some())
        .filter_map(|(atom, _)| u32::try_from(atom).ok())
        .collect();
    let count = matched.len();
    count_kept(
        reference,
        &AtomSelection::from_sorted(matched),
        count,
        |atom| target_of.get(atom as usize).copied().flatten(),
        target_positions,
        (cutoff, tolerance),
        backend,
        context,
    )
}

fn count_kept(
    reference: &Structure,
    selection: &AtomSelection,
    matched_atoms: usize,
    target_of: impl Fn(u32) -> Option<u32> + Sync,
    target_positions: &[[f32; 3]],
    (cutoff, tolerance): (f32, f32),
    backend: SpatialBackend,
    context: &ExecutionContext,
) -> Result<NativeContacts, NativeError> {
    let kept_cutoff = f64::from(cutoff) + f64::from(tolerance);
    let kept_cutoff_squared = kept_cutoff * kept_cutoff * (1.0 + DISTANCE_ROUNDING);

    // Both totals are sums over pairs, so each block keeps its own counters and
    // nothing pair-proportional is ever retained. Integer addition is exact and
    // the blocks merge in plan order, so the result is worker-count independent.
    let query = PairQuery {
        positions: reference.positions(),
        left: selection,
        right: selection,
        cutoff,
        options: SpatialSearchOptions::with_backend(backend),
        periodic: None,
        context,
    };
    let counts = reduce_pairs_within_unsorted(
        &query,
        || (0usize, 0usize),
        |totals, pair| {
            totals.0 += 1;
            let (Some(first), Some(second)) = (target_of(pair.first), target_of(pair.second))
            else {
                return;
            };
            let (Some(&a), Some(&b)) = (
                target_positions.get(first as usize),
                target_positions.get(second as usize),
            ) else {
                return;
            };
            if molframe_geom::distance_squared(a, b) <= kept_cutoff_squared {
                totals.1 += 1;
            }
        },
    )?;

    let mut native = 0usize;
    let mut kept = 0usize;
    for (block_native, block_kept) in counts {
        native += block_native;
        kept += block_kept;
    }

    let fraction = if native == 0 {
        1.0
    } else {
        usize_to_f64(kept) / usize_to_f64(native)
    };
    Ok(NativeContacts {
        native,
        kept,
        matched_atoms,
        fraction,
    })
}

#[cfg(test)]
#[path = "native_tests.rs"]
mod tests;

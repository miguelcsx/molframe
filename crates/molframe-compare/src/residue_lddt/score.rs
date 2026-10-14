//! The residue-aware lDDT over a table of reference atoms.
//!
//! The reference defines the tested pairs: every pair of its heavy atoms from
//! different residues within the inclusion radius. A pair is scored by the
//! fraction of tolerances its model distance meets, and counts as failed when
//! either atom is absent from the model. Every pair is credited to both of its
//! residues, and the global score is the same pooled fraction over all pairs, so
//! a residue's score and the global score weigh a distance identically. A
//! fixed-radius cell list visits only candidate pairs, `O(n + k)` for `k` pairs.

use super::options::{ResidueLddt, ResidueLddtOptions, ResidueScore, SymmetryPolicy};
use super::resolve::{Outcome, resolve};
use super::table::{ResidueCorrespondence, Table};
use crate::local_distance::conservative_f32;
use crate::mapped::MappedCompareError;
use crate::numeric::{u64_to_f64, usize_to_f64};
use crate::{CompareError, EmptyLddtPolicy};
use molframe_chem::ComponentProvider;
use molframe_core::{AtomSelection, ExecutionContext};
use molframe_spatial::{
    PairQuery, SpatialBackend, SpatialError, SpatialSearchOptions, reduce_pairs_within_unsorted,
};
use std::collections::BTreeMap;

/// Scores `model` against `reference` over the residues in `residues`.
///
/// The reference residues named in `residues` form the domain; atoms of other
/// residues, hydrogens and the reference's own waters never enter. Component
/// definitions come from `provider` and only decide which atoms are
/// interchangeable.
///
/// # Errors
///
/// Returns a mapping error for a residue index outside either structure or a
/// provider failure, and a comparison error for invalid options, non-finite
/// coordinates, or an empty domain under [`EmptyLddtPolicy::Error`].
pub fn residue_lddt(
    reference: &molframe_core::structure::Structure,
    model: &molframe_core::structure::Structure,
    residues: &[ResidueCorrespondence],
    provider: &dyn ComponentProvider,
    options: &ResidueLddtOptions,
    context: &ExecutionContext,
) -> Result<ResidueLddt, MappedCompareError> {
    options.validate()?;
    let table = Table::build(reference, model, residues)?;
    let Outcome {
        assigned,
        renamed,
        unresolved,
    } = match options.symmetry {
        SymmetryPolicy::Disabled => Outcome {
            assigned: table.identity(),
            renamed: 0,
            unresolved: 0,
        },
        SymmetryPolicy::Resolve { automorphism_limit } => {
            resolve(&table, provider, options, automorphism_limit, context)?
        }
    };
    let tallies = tally(&table, &assigned, options, context)?;

    let mut considered = 0;
    let mut hits = 0;
    let mut scores = Vec::with_capacity(residues.len());
    for (pair, (seen, passed)) in residues.iter().zip(&tallies) {
        considered += seen;
        hits += passed;
        scores.push(ResidueScore {
            reference: pair.reference,
            model: pair.model,
            considered: *seen,
            preserved: *passed,
            score: fraction(*passed, *seen, options),
        });
    }
    // Each pair was credited to both residues.
    let pairs = considered / 2;
    let score = match fraction(hits, considered, options) {
        Some(score) => score,
        None => match options.empty_policy {
            EmptyLddtPolicy::Perfect => 1.0,
            EmptyLddtPolicy::Error => return Err(CompareError::NoComparablePairs.into()),
        },
    };
    Ok(ResidueLddt {
        score,
        considered: pairs,
        residues: scores,
        renamed_residues: renamed,
        unresolved_symmetry: unresolved,
    })
}

fn fraction(hits: u64, considered: u64, options: &ResidueLddtOptions) -> Option<f64> {
    if considered == 0 {
        return None;
    }
    Some(u64_to_f64(hits) / (u64_to_f64(considered) * usize_to_f64(options.tolerances.len())))
}

/// How many tolerances a model distance meets for a reference distance.
pub(super) fn preserved(
    tolerances: &[f64],
    reference_distance: f64,
    first: [f32; 3],
    second: [f32; 3],
) -> u64 {
    let error = (molframe_geom::distance(first, second) - reference_distance).abs();
    tolerances
        .iter()
        .filter(|tolerance| error < **tolerance)
        .count() as u64
}

/// Reference pairs within the radius that pass the residue rules and `keep`.
pub(super) fn eligible_pairs(
    table: &Table<'_>,
    options: &ResidueLddtOptions,
    context: &ExecutionContext,
    keep: impl Fn(u32, u32) -> bool + Sync,
) -> Result<Vec<(u32, u32, f32)>, MappedCompareError> {
    let blocks = fold_pairs(
        table,
        options,
        context,
        Vec::new,
        |rows, first, second, d| {
            if keep(first, second) {
                // Stored narrow on purpose: rows only feed per-residue tallies.
                #[allow(clippy::cast_possible_truncation)]
                rows.push((first, second, d as f32));
            }
        },
    )?;
    Ok(blocks.into_iter().flatten().collect())
}

fn tally(
    table: &Table<'_>,
    assigned: &[u32],
    options: &ResidueLddtOptions,
    context: &ExecutionContext,
) -> Result<Vec<(u64, u64)>, MappedCompareError> {
    let blocks = fold_pairs(
        table,
        options,
        context,
        BTreeMap::<u32, (u64, u64)>::new,
        |counts, first, second, distance| {
            let hits = match (
                table.model[assigned[first as usize] as usize],
                table.model[assigned[second as usize] as usize],
            ) {
                (Some(left), Some(right)) => preserved(&options.tolerances, distance, left, right),
                _ => 0,
            };
            for owner in [table.owner[first as usize], table.owner[second as usize]] {
                let entry = counts.entry(owner).or_insert((0, 0));
                entry.0 += 1;
                entry.1 += hits;
            }
        },
    )?;
    let mut totals = vec![(0, 0); table.spans.len()];
    for block in blocks {
        for (owner, (seen, passed)) in block {
            totals[owner as usize].0 += seen;
            totals[owner as usize].1 += passed;
        }
    }
    Ok(totals)
}

/// Folds every testable reference pair into per-block state, in block order.
fn fold_pairs<T: Send>(
    table: &Table<'_>,
    options: &ResidueLddtOptions,
    context: &ExecutionContext,
    init: impl Fn() -> T + Sync,
    fold: impl Fn(&mut T, u32, u32, f64) + Sync,
) -> Result<Vec<T>, MappedCompareError> {
    let count = u32::try_from(table.reference.len())
        .map_err(|_| CompareError::Spatial(SpatialError::NumericRangeExceeded))?;
    let cutoff =
        conservative_f32(options.inclusion_radius).ok_or(CompareError::InvalidDistanceCutoff)?;
    let selection = AtomSelection::All(count);
    let limit = options.inclusion_radius * options.inclusion_radius;
    let query = PairQuery {
        positions: &table.reference,
        left: &selection,
        right: &selection,
        cutoff,
        options: SpatialSearchOptions::with_backend(SpatialBackend::CellList),
        periodic: None,
        context,
    };
    let blocks = reduce_pairs_within_unsorted(&query, init, |state, pair| {
        let (first, second) = (pair.first, pair.second);
        if !table.testable(
            table.owner[first as usize],
            table.owner[second as usize],
            options.sequence_separation,
        ) {
            return;
        }
        let distance_squared = molframe_geom::distance_squared(
            table.reference[first as usize],
            table.reference[second as usize],
        );
        if distance_squared <= limit {
            fold(state, first, second, distance_squared.sqrt());
        }
    })
    .map_err(CompareError::from)?;
    Ok(blocks)
}

#[cfg(test)]
#[path = "score_tests.rs"]
mod tests;

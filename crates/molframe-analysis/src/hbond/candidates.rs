//! Conservative oriented hydrogen-bond candidates for linear coordinate motion.
//!
//! A spatial query uses a displacement envelope after removing a common
//! translation. Preparation costs O(atoms + spatial candidates); no coordinates
//! are copied. Exact distance and angle tests belong to the sampled pose.

use super::{
    HydrogenBondError, donor_hydrogens, eligible_hydrogen, eligible_pair, orientations,
    role_selection,
};
use molframe_core::selection::AtomSelection;
use molframe_core::{AtomIndex, ExecutionContext, Structure};
use molframe_spatial::{
    PairQuery, SpatialBackend, SpatialSearchOptions, reduce_pairs_within_unsorted,
};

/// A chemically eligible triple; it is not yet a geometrically accepted bond.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HydrogenBondCandidate {
    /// CCD-perceived donor atom.
    pub donor: AtomIndex,
    /// Explicit hydrogen covalently attached to the donor.
    pub hydrogen: AtomIndex,
    /// CCD-perceived acceptor atom.
    pub acceptor: AtomIndex,
}

define_soa_table! {
    /// Topology-backed oriented rows that may satisfy geometry during a segment.
    pub struct HydrogenBondCandidateTable for HydrogenBondCandidate {
        /// Donor atom indices.
        donor: AtomIndex,
        /// Explicit hydrogen atom indices.
        hydrogen: AtomIndex,
        /// Acceptor atom indices.
        acceptor: AtomIndex,
    }
}

/// Nonperiodic motion envelope and an explicit bound on candidate storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HydrogenBondCandidateOptions {
    /// Maximum accepted donor–acceptor distance in ångström.
    pub maximum_donor_acceptor_distance: f32,
    /// Maximum oriented triples; exceeding it rejects the whole result.
    pub maximum_candidates: usize,
    /// Spatial implementation used for the conservative pair query.
    pub backend: SpatialBackend,
}

impl HydrogenBondCandidateOptions {
    /// Validates the physical distance envelope and addressable storage limit.
    ///
    /// # Errors
    /// Rejects nonpositive or nonfinite distances and budgets outside one to u32.
    pub fn validate(self) -> Result<(), HydrogenBondError> {
        if !self.maximum_donor_acceptor_distance.is_finite()
            || self.maximum_donor_acceptor_distance <= 0.0
            || self.maximum_candidates == 0
            || self.maximum_candidates > u32::MAX as usize
        {
            return Err(HydrogenBondError::InvalidOptions);
        }
        Ok(())
    }
}

#[derive(Default)]
struct Candidates {
    rows: Vec<HydrogenBondCandidate>,
    exceeded: bool,
}

/// Finds every eligible triple whose heavy atoms may approach within the
/// distance cutoff during linear interpolation between two borrowed columns.
///
/// The output is conservative: angle and instantaneous distance are deliberately
/// not filtered here, so contacts can appear and disappear without rebuilding
/// topology. Roles, covalent attachment and alternate-location compatibility
/// use the same chemistry predicates as static hydrogen-bond detection.
/// Results are sorted by `(donor, hydrogen, acceptor)` and contain no duplicates.
/// Coordinates and distances are Cartesian, without periodic images.
///
/// # Errors
/// Rejects invalid cutoffs, budgets, nonfinite or incorrectly sized coordinate
/// columns, absent chemistry, spatial failures and excessive candidate counts.
pub fn hydrogen_bond_candidates_over_segment(
    structure: &Structure,
    start: &[[f32; 3]],
    end: &[[f32; 3]],
    options: HydrogenBondCandidateOptions,
    context: &ExecutionContext,
) -> Result<HydrogenBondCandidateTable, HydrogenBondError> {
    validate(structure, start, end, options)?;
    if !structure.data().bonds.is_available() {
        return Err(HydrogenBondError::MissingChemistry);
    }
    let donors = role_selection(structure, molframe_core::HBOND_DONOR_ANNOTATION)?;
    let acceptors = role_selection(structure, molframe_core::HBOND_ACCEPTOR_ANNOTATION)?;
    let cutoff = envelope_cutoff(
        start,
        end,
        &donors,
        &acceptors,
        options.maximum_donor_acceptor_distance,
    )?;
    let adjacency = structure.data().bonds.adjacency(structure.atom_count());
    let query = PairQuery {
        positions: start,
        left: &donors,
        right: &acceptors,
        cutoff,
        options: SpatialSearchOptions::with_backend(options.backend),
        periodic: None,
        context,
    };
    let parts =
        reduce_pairs_within_unsorted(&query, Candidates::default, |out: &mut Candidates, pair| {
            if out.exceeded {
                return;
            }
            for (donor, acceptor) in orientations(pair.first, pair.second, &donors, &acceptors) {
                if !eligible_pair(structure, adjacency, donor, acceptor) {
                    continue;
                }
                if !approaches(
                    start,
                    end,
                    donor,
                    acceptor,
                    options.maximum_donor_acceptor_distance,
                ) {
                    continue;
                }
                for hydrogen in donor_hydrogens(structure, adjacency, donor) {
                    if !eligible_hydrogen(structure, donor, hydrogen, acceptor) {
                        continue;
                    }
                    if out.rows.len() == options.maximum_candidates {
                        out.exceeded = true;
                        return;
                    }
                    out.rows.push(HydrogenBondCandidate {
                        donor: AtomIndex::new(donor),
                        hydrogen: AtomIndex::new(hydrogen),
                        acceptor: AtomIndex::new(acceptor),
                    });
                }
            }
        })?;
    let mut rows = Vec::new();
    for part in parts {
        if part.exceeded || part.rows.len() > options.maximum_candidates - rows.len() {
            return Err(HydrogenBondError::CandidateBudget {
                maximum: options.maximum_candidates,
            });
        }
        rows.extend(part.rows);
    }
    rows.sort_unstable_by_key(|row| (row.donor.get(), row.hydrogen.get(), row.acceptor.get()));
    rows.dedup();
    Ok(rows.into_iter().collect())
}

fn validate(
    structure: &Structure,
    start: &[[f32; 3]],
    end: &[[f32; 3]],
    options: HydrogenBondCandidateOptions,
) -> Result<(), HydrogenBondError> {
    options.validate()?;
    if start.len() != structure.atom_count() as usize
        || end.len() != start.len()
        || start
            .iter()
            .chain(end)
            .flatten()
            .any(|value| !value.is_finite())
    {
        return Err(HydrogenBondError::InvalidOptions);
    }
    Ok(())
}

fn envelope_cutoff(
    start: &[[f32; 3]],
    end: &[[f32; 3]],
    donors: &AtomSelection,
    acceptors: &AtomSelection,
    distance: f32,
) -> Result<f32, HydrogenBondError> {
    let drift = match donors.iter().chain(acceptors.iter()).next() {
        Some(row) => std::array::from_fn(|i| {
            f64::from(end[row as usize][i]) - f64::from(start[row as usize][i])
        }),
        _ => [0.0; 3],
    };
    let mut maximum_squared = 0.0f64;
    for row in donors.iter().chain(acceptors.iter()) {
        let first = start[row as usize];
        let last = end[row as usize];
        let displacement: [f64; 3] =
            std::array::from_fn(|i| f64::from(last[i]) - f64::from(first[i]) - drift[i]);
        let squared = displacement.into_iter().map(|value| value * value).sum();
        maximum_squared = maximum_squared.max(squared);
    }
    let cutoff = f64::from(distance) + 2.0 * maximum_squared.sqrt();
    if cutoff >= f64::from(f32::MAX) {
        return Err(HydrogenBondError::InvalidOptions);
    }
    Ok(crate::numeric::f64_to_f32(cutoff).next_up())
}

fn approaches(
    start: &[[f32; 3]],
    end: &[[f32; 3]],
    donor: u32,
    acceptor: u32,
    cutoff: f32,
) -> bool {
    let relative = |positions: &[[f32; 3]]| -> [f64; 3] {
        std::array::from_fn(|i| {
            f64::from(positions[donor as usize][i]) - f64::from(positions[acceptor as usize][i])
        })
    };
    let first = relative(start);
    let last = relative(end);
    let velocity: [f64; 3] = std::array::from_fn(|i| last[i] - first[i]);
    let squared_speed: f64 = velocity.iter().map(|value| value * value).sum();
    let projection: f64 = first
        .iter()
        .zip(velocity)
        .map(|(position, speed)| position * speed)
        .sum();
    let amount = if squared_speed > 0.0 {
        (-projection / squared_speed).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let minimum_squared: f64 = first
        .iter()
        .zip(velocity)
        .map(|(position, speed)| (position + speed * amount).powi(2))
        .sum();
    let limit = f64::from(cutoff).powi(2);
    // Include floating-point boundary contacts in the candidate superset.
    minimum_squared <= limit + limit * 1e-6
}

#[cfg(test)]
#[path = "candidates_tests.rs"]
mod tests;

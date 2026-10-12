//! Sampled hydrogen-bond snapshots cost O(atoms + candidates) without copying coordinates.

use super::{
    HydrogenBond, HydrogenBondCandidateTable, HydrogenBondError, HydrogenBondOptions,
    HydrogenBondPolicy, HydrogenBondTable, eligible_hydrogen, eligible_pair, measure_positions,
    role_selection, validate_options,
};
use molframe_core::{CoordinateSample, ExecutionContext, Structure};
use molframe_spatial::SpatialError;

/// Applies explicit-hydrogen geometry to a chemically eligible candidate table
/// at one borrowed Cartesian sample. Output is sorted by donor, hydrogen and
/// acceptor, independent of candidate ordering, and contains no duplicates.
///
/// The table must cover this sample's resident interval. Preparation and its
/// interval identity and a distance envelope at least as large as the requested
/// cutoff remain the caller's responsibility; every supplied triple
/// is checked against the original molecule's chemical roles and connectivity.
/// Periodic images are refused because candidates describe nonperiodic motion.
///
/// # Errors
/// Refuses invalid geometry, periodic requests, absent chemistry, mismatched or
/// nonfinite coordinates, invalid chemical triples, cancellation and memory exhaustion.
pub fn hydrogen_bonds_at_sample(
    structure: &Structure,
    sample: CoordinateSample<'_>,
    candidates: &HydrogenBondCandidateTable,
    options: HydrogenBondOptions,
    policy: HydrogenBondPolicy,
    context: &ExecutionContext,
) -> Result<HydrogenBondTable, HydrogenBondError> {
    validate_options(options, policy)?;
    if options.periodic {
        return Err(HydrogenBondError::InvalidOptions);
    }
    if sample.len() != structure.atom_count() as usize {
        return Err(HydrogenBondError::InvalidSample);
    }
    if !structure.data().bonds.is_available() {
        return Err(HydrogenBondError::MissingChemistry);
    }
    check(context)?;
    let donors = role_selection(structure, molframe_core::HBOND_DONOR_ANNOTATION)?;
    let acceptors = role_selection(structure, molframe_core::HBOND_ACCEPTOR_ANNOTATION)?;
    for row in 0..sample.len() {
        check(context)?;
        if sample.position(row).is_none() {
            return Err(HydrogenBondError::InvalidSample);
        }
    }
    let bytes = candidates
        .len()
        .checked_mul(2 * std::mem::size_of::<HydrogenBond>())
        .ok_or(HydrogenBondError::InvalidSample)?;
    let _memory = context.try_reserve(bytes).map_err(SpatialError::from)?;
    let adjacency = structure.data().bonds.adjacency(structure.atom_count());
    let mut output = Vec::with_capacity(candidates.len());
    for row in candidates.iter() {
        check(context)?;
        let donor = row.donor.get();
        let hydrogen = row.hydrogen.get();
        let acceptor = row.acceptor.get();
        if !donors.contains(donor)
            || !acceptors.contains(acceptor)
            || !eligible_pair(structure, adjacency, donor, acceptor)
            || adjacency
                .neighbours(row.donor)
                .binary_search(&row.hydrogen)
                .is_err()
            || !structure
                .data()
                .atom(row.hydrogen)
                .and_then(molframe_core::structure::AtomRef::element)
                .is_some_and(molframe_core::Element::is_hydrogen)
            || !eligible_hydrogen(structure, donor, hydrogen, acceptor)
        {
            return Err(HydrogenBondError::InvalidSample);
        }
        let position = |atom: u32| {
            sample
                .position(atom as usize)
                .ok_or(HydrogenBondError::InvalidSample)
        };
        let Some(bond) = measure_positions(
            [donor, hydrogen, acceptor],
            [position(donor)?, position(hydrogen)?, position(acceptor)?],
            None,
        ) else {
            continue;
        };
        if bond.donor_acceptor_distance <= options.maximum_donor_acceptor_distance
            && bond.angle_degrees >= options.minimum_angle_degrees
            && policy
                .maximum_hydrogen_acceptor_distance
                .is_none_or(|limit| f64::from(bond.hydrogen_acceptor_distance) <= limit)
        {
            output.push(bond);
        }
    }
    output.sort_unstable_by_key(|bond| (bond.donor, bond.hydrogen, bond.acceptor));
    output.dedup_by_key(|bond| (bond.donor, bond.hydrogen, bond.acceptor));
    Ok(output.into_iter().collect())
}

fn check(context: &ExecutionContext) -> Result<(), HydrogenBondError> {
    context
        .cancellation()
        .check()
        .map_err(|_| SpatialError::Cancelled)?;
    Ok(())
}

#[cfg(test)]
#[path = "sample_tests.rs"]
mod tests;

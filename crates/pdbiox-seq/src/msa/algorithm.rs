//! Guide-tree construction, progressive merge and deterministic refinement.

use super::profile::{Profile, align};
use super::{MsaError, MsaOptions};
use crate::Scoring;
use crate::align::global;

/// Progressively aligns sequences from a deterministic average-linkage guide tree.
///
/// Pairwise affine alignments define normalized mismatch distances. The closest
/// two clusters are repeatedly merged, and each guide-tree join performs an
/// affine profile-profile alignment. Optional leave-one-out passes keep a
/// realignment only when it improves the full sum-of-pairs objective.
///
/// # Errors
///
/// Returns [`MsaError`] for gap scores that reward gaps or dimensions that
/// cannot be represented safely.
pub fn progressive_msa(sequences: &[&[u8]], options: MsaOptions) -> Result<Vec<Vec<u8>>, MsaError> {
    validate(sequences, options)?;
    if sequences.is_empty() {
        return Ok(Vec::new());
    }
    let distances = pairwise_distances(sequences, options.scoring)?;
    let mut profiles = sequences
        .iter()
        .enumerate()
        .map(|(index, sequence)| Some(Profile::singleton(index, sequence)))
        .collect::<Vec<_>>();
    let mut members = (0..sequences.len())
        .map(|index| Some(vec![index]))
        .collect::<Vec<_>>();
    while active_count(&profiles) > 1 {
        let Some((left, right)) = closest_pair(&members, &distances)? else {
            return Err(MsaError::DimensionOverflow);
        };
        let Some(left_profile) = profiles[left].take() else {
            return Err(MsaError::DimensionOverflow);
        };
        let Some(right_profile) = profiles[right].take() else {
            return Err(MsaError::DimensionOverflow);
        };
        profiles[left] = Some(align(&left_profile, &right_profile, options.scoring)?);
        let Some(mut right_members) = members[right].take() else {
            return Err(MsaError::DimensionOverflow);
        };
        let Some(left_members) = members[left].as_mut() else {
            return Err(MsaError::DimensionOverflow);
        };
        left_members.append(&mut right_members);
        left_members.sort_unstable();
    }
    let Some(mut profile) = profiles.into_iter().flatten().next() else {
        return Err(MsaError::DimensionOverflow);
    };
    for _ in 0..options.refinement_passes {
        profile = refine(profile, sequences, options.scoring)?;
    }
    profile.rows.sort_by_key(|(index, _)| *index);
    Ok(profile.rows.into_iter().map(|(_, row)| row).collect())
}

/// Convenience form of [`progressive_msa`] without iterative refinement.
///
/// The caller still supplies the complete scoring policy; only the optional
/// refinement count is selected by this convenience entry point.
///
/// # Errors
///
/// Returns [`MsaError`] under the same conditions as [`progressive_msa`].
pub fn msa(sequences: &[&[u8]], scoring: Scoring) -> Result<Vec<Vec<u8>>, MsaError> {
    progressive_msa(sequences, MsaOptions::progressive(scoring))
}

fn validate(sequences: &[&[u8]], options: MsaOptions) -> Result<(), MsaError> {
    if options.scoring.gap_open > 0 || options.scoring.gap_extend > 0 {
        return Err(MsaError::InvalidGapScore);
    }
    let total = sequences
        .iter()
        .try_fold(0usize, |sum, sequence| sum.checked_add(sequence.len()));
    let Some(total) = total else {
        return Err(MsaError::DimensionOverflow);
    };
    let matrix_fits = total
        .checked_add(1)
        .and_then(|side| side.checked_mul(side))
        .is_some();
    let counts_fit = u32::try_from(total).is_ok()
        && u32::try_from(sequences.len()).is_ok()
        && sequences
            .iter()
            .all(|sequence| u32::try_from(sequence.len()).is_ok());
    if matrix_fits && counts_fit {
        Ok(())
    } else {
        Err(MsaError::DimensionOverflow)
    }
}

fn pairwise_distances(sequences: &[&[u8]], scoring: Scoring) -> Result<Vec<Vec<f64>>, MsaError> {
    let mut distances = vec![vec![0.0; sequences.len()]; sequences.len()];
    for left in 0..sequences.len() {
        for right in left + 1..sequences.len() {
            let alignment = global(sequences[left], sequences[right], scoring)?;
            let compared = alignment
                .columns
                .iter()
                .filter(|column| column.left.is_some() && column.right.is_some())
                .count();
            let matches = alignment
                .columns
                .iter()
                .filter(|column| match (column.left, column.right) {
                    (Some(first), Some(second)) => {
                        sequences[left][first] == sequences[right][second]
                    }
                    _ => false,
                })
                .count();
            let normalizer = alignment.columns.len().max(1);
            let normalizer = count_as_f64(normalizer)?;
            let identity = count_as_f64(matches)? / normalizer;
            let coverage = count_as_f64(compared)? / normalizer;
            let distance = 1.0 - identity * coverage;
            distances[left][right] = distance;
            distances[right][left] = distance;
        }
    }
    Ok(distances)
}

fn active_count(profiles: &[Option<Profile>]) -> usize {
    profiles.iter().filter(|profile| profile.is_some()).count()
}

fn closest_pair(
    members: &[Option<Vec<usize>>],
    distances: &[Vec<f64>],
) -> Result<Option<(usize, usize)>, MsaError> {
    let mut best: Option<(f64, usize, usize)> = None;
    for left in 0..members.len() {
        let Some(left_members) = &members[left] else {
            continue;
        };
        for (right, right_members) in members.iter().enumerate().skip(left + 1) {
            let Some(right_members) = right_members else {
                continue;
            };
            let distance = average_distance(left_members, right_members, distances)?;
            match best {
                Some((current, _, _)) if distance >= current => {}
                _ => best = Some((distance, left, right)),
            }
        }
    }
    Ok(best.map(|(_, left, right)| (left, right)))
}

fn average_distance(
    left: &[usize],
    right: &[usize],
    distances: &[Vec<f64>],
) -> Result<f64, MsaError> {
    let mut total = 0.0;
    for &first in left {
        for &second in right {
            total += distances[first][second];
        }
    }
    let pairs = left
        .len()
        .checked_mul(right.len())
        .ok_or(MsaError::DimensionOverflow)?;
    Ok(total / count_as_f64(pairs)?)
}

fn count_as_f64(value: usize) -> Result<f64, MsaError> {
    u32::try_from(value)
        .map(f64::from)
        .map_err(|_| MsaError::DimensionOverflow)
}

fn refine(
    mut profile: Profile,
    sequences: &[&[u8]],
    scoring: Scoring,
) -> Result<Profile, MsaError> {
    for (index, sequence) in sequences.iter().enumerate() {
        if profile.rows.len() <= 1 {
            break;
        }
        let previous_score = profile.score(scoring)?;
        let remainder = profile.without_row(index);
        let candidate = align(&remainder, &Profile::singleton(index, sequence), scoring)?;
        if candidate.score(scoring)? > previous_score {
            profile = candidate;
        }
    }
    Ok(profile)
}

#[cfg(test)]
#[path = "algorithm_tests.rs"]
mod tests;

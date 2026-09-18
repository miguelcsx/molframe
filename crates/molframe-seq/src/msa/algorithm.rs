//! Guide-tree construction, progressive merge and deterministic refinement.

use super::guide::single_linkage_guide;
use super::profile::{Profile, align, alignment_peak_bytes};
use super::{MsaError, MsaOptions};
use crate::Scoring;

/// Progressively aligns sequences from a deterministic single-linkage guide tree.
///
/// Bit-parallel longest-common-subsequence scores define normalized guide
/// distances. An incremental minimum spanning tree keeps guide storage linear
/// in the sequence count, and each tree join performs an affine
/// profile-profile alignment. Optional leave-one-out passes keep a realignment
/// only when it improves the full sum-of-pairs objective.
///
/// # Errors
///
/// Returns [`MsaError`] for gap scores that reward gaps, dimensions that cannot
/// be represented safely, or work that exceeds the explicit memory ceiling.
pub fn progressive_msa(sequences: &[&[u8]], options: MsaOptions) -> Result<Vec<Vec<u8>>, MsaError> {
    validate(sequences, options)?;
    if sequences.is_empty() {
        return Ok(Vec::new());
    }
    let guide = single_linkage_guide(sequences, options.memory_limit_bytes)?;
    let guide_storage = guide
        .capacity()
        .checked_mul(size_of::<(usize, usize)>())
        .ok_or(MsaError::DimensionOverflow)?;
    let mut profiles = sequences
        .iter()
        .enumerate()
        .map(|(index, sequence)| Some(Profile::singleton(index, sequence)))
        .collect::<Vec<_>>();
    check_memory(
        profile_storage(&profiles)?
            .checked_add(guide_storage)
            .ok_or(MsaError::DimensionOverflow)?,
        options.memory_limit_bytes,
    )?;
    for (left, right) in guide {
        let Some(left_profile) = profiles[left].take() else {
            return Err(MsaError::DimensionOverflow);
        };
        let Some(right_profile) = profiles[right].take() else {
            return Err(MsaError::DimensionOverflow);
        };
        let alignment_bytes = alignment_peak_bytes(&left_profile, &right_profile)?;
        let required = profile_storage(&profiles)?
            .checked_add(guide_storage)
            .and_then(|bytes| bytes.checked_add(alignment_bytes))
            .ok_or(MsaError::DimensionOverflow)?;
        check_memory(required, options.memory_limit_bytes)?;
        profiles[left] = Some(align(
            &left_profile,
            &right_profile,
            options.scoring,
            options.memory_limit_bytes,
        )?);
    }
    let Some(mut profile) = profiles.into_iter().flatten().next() else {
        return Err(MsaError::DimensionOverflow);
    };
    for _ in 0..options.refinement_passes {
        profile = refine(profile, sequences, options)?;
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
    if options.memory_limit_bytes == 0 {
        return Err(MsaError::InvalidMemoryLimit {
            requested: options.memory_limit_bytes,
        });
    }
    let total = sequences
        .iter()
        .try_fold(0usize, |sum, sequence| sum.checked_add(sequence.len()));
    let Some(total) = total else {
        return Err(MsaError::DimensionOverflow);
    };
    let counts_fit = u32::try_from(total).is_ok()
        && u32::try_from(sequences.len()).is_ok()
        && sequences
            .iter()
            .all(|sequence| u32::try_from(sequence.len()).is_ok());
    if counts_fit {
        Ok(())
    } else {
        Err(MsaError::DimensionOverflow)
    }
}

fn profile_storage(profiles: &[Option<Profile>]) -> Result<usize, MsaError> {
    profiles
        .iter()
        .flatten()
        .try_fold(0_usize, |total, profile| {
            profile
                .storage_bytes()
                .and_then(|bytes| total.checked_add(bytes))
        })
        .ok_or(MsaError::DimensionOverflow)
}

fn check_memory(required: usize, limit: usize) -> Result<(), MsaError> {
    if required <= limit {
        Ok(())
    } else {
        Err(MsaError::MemoryLimit { required, limit })
    }
}

fn refine(
    mut profile: Profile,
    sequences: &[&[u8]],
    options: MsaOptions,
) -> Result<Profile, MsaError> {
    let mut current_score = profile.score(options.scoring)?;
    for (index, sequence) in sequences.iter().enumerate() {
        if profile.rows.len() <= 1 {
            break;
        }
        let profile_bytes = profile.storage_bytes().ok_or(MsaError::DimensionOverflow)?;
        let cloning_peak = profile_bytes
            .checked_mul(2)
            .and_then(|bytes| bytes.checked_add(sequence.len()))
            .ok_or(MsaError::DimensionOverflow)?;
        check_memory(cloning_peak, options.memory_limit_bytes)?;
        let remainder = profile.without_row(index);
        let singleton = Profile::singleton(index, sequence);
        let required = profile_bytes
            .checked_add(alignment_peak_bytes(&remainder, &singleton)?)
            .ok_or(MsaError::DimensionOverflow)?;
        check_memory(required, options.memory_limit_bytes)?;
        let candidate = align(
            &remainder,
            &singleton,
            options.scoring,
            options.memory_limit_bytes,
        )?;
        let candidate_score = candidate.score(options.scoring)?;
        if candidate_score > current_score {
            profile = candidate;
            current_score = candidate_score;
        }
    }
    Ok(profile)
}

#[cfg(test)]
#[path = "algorithm_tests.rs"]
mod tests;

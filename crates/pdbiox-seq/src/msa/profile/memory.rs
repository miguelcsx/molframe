//! Conservative allocation accounting for profile dynamic programming.

use super::Profile;
use crate::MsaError;

pub(in crate::msa) fn alignment_peak_bytes(
    left: &Profile,
    right: &Profile,
) -> Result<usize, MsaError> {
    let left_width = left.width();
    let right_width = right.width();
    let trace_cells = left_width
        .checked_mul(right_width)
        .ok_or(MsaError::DimensionOverflow)?;
    let trace = trace_cells
        .checked_add(1)
        .map(|cells| cells / 2)
        .ok_or(MsaError::DimensionOverflow)?;
    let score_rows = right_width
        .checked_add(1)
        .and_then(|width| width.checked_mul(6 * size_of::<i64>()))
        .ok_or(MsaError::DimensionOverflow)?;
    let column_stats = stats_bytes(left)?
        .checked_add(stats_bytes(right)?)
        .ok_or(MsaError::DimensionOverflow)?;
    let gap_scores = left_width
        .checked_add(right_width)
        .and_then(|width| width.checked_mul(2 * size_of::<i64>()))
        .ok_or(MsaError::DimensionOverflow)?;
    let projected_width = left_width
        .checked_add(right_width)
        .ok_or(MsaError::DimensionOverflow)?;
    let traceback_columns = projected_width
        .checked_mul(2 * size_of::<Option<usize>>())
        .ok_or(MsaError::DimensionOverflow)?;
    let output_rows = left
        .rows
        .len()
        .checked_add(right.rows.len())
        .ok_or(MsaError::DimensionOverflow)?;
    let output = output_rows
        .checked_mul(projected_width)
        .and_then(|bytes| {
            output_rows
                .checked_mul(size_of::<(usize, Vec<u8>)>())
                .and_then(|headers| bytes.checked_add(headers))
        })
        .ok_or(MsaError::DimensionOverflow)?;
    let inputs = left
        .storage_bytes()
        .and_then(|bytes| right.storage_bytes()?.checked_add(bytes))
        .ok_or(MsaError::DimensionOverflow)?;
    [
        inputs,
        trace,
        score_rows,
        column_stats,
        gap_scores,
        traceback_columns,
        output,
    ]
    .into_iter()
    .try_fold(0_usize, usize::checked_add)
    .ok_or(MsaError::DimensionOverflow)
}

fn stats_bytes(profile: &Profile) -> Result<usize, MsaError> {
    let width = profile.width();
    let cells = width
        .checked_mul(profile.rows.len())
        .ok_or(MsaError::DimensionOverflow)?;
    let sparse = cells
        .checked_mul(size_of::<u8>() + size_of::<u32>())
        .ok_or(MsaError::DimensionOverflow)?;
    width
        .checked_add(1)
        .and_then(|offsets| offsets.checked_mul(size_of::<u32>()))
        .and_then(|offsets| {
            width
                .checked_mul(size_of::<u32>())
                .and_then(|residues| offsets.checked_add(residues))
        })
        .and_then(|fixed| fixed.checked_add(sparse))
        .ok_or(MsaError::DimensionOverflow)
}

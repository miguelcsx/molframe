//! Working-set and identity checks for bounded CIF batches.

use pdbiox_core::{BatchDemand, Diagnostic, StructureBatchError};

const RETAINED_BYTES_PER_ROW: usize = 128;
const DICTIONARY_HEADROOM: usize = 64 * 1024;

pub(super) fn row_capacity(demand: BatchDemand) -> Result<u32, StructureBatchError> {
    let rows = demand
        .max_bytes
        .saturating_sub(DICTIONARY_HEADROOM)
        .checked_div(RETAINED_BYTES_PER_ROW)
        .map_or(0, |value| value)
        .min(demand.max_rows)
        .min(u32::MAX as usize);
    if rows == 0 {
        return Err(StructureBatchError::DemandTooSmall {
            required: DICTIONARY_HEADROOM + RETAINED_BYTES_PER_ROW,
            available: demand.max_bytes,
        });
    }
    u32::try_from(rows).map_err(|_| identity_overflow())
}

pub(super) const fn reserve_bytes(_capacity: u32, demand: BatchDemand) -> usize {
    demand.max_bytes
}

pub(super) fn identity_overflow() -> StructureBatchError {
    StructureBatchError::Diagnostic(
        Diagnostic::new(pdbiox_core::Code::E1903).with_message("CIF batch identity exceeds u64"),
    )
}

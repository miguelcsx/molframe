//! Overflow-safe working-set estimates for incremental BCIF batches.

use super::{DICTIONARY_HEADROOM, RETAINED_BYTES_PER_ROW};
use pdbiox_core::{Code, Diagnostic, StructureBatchError};

pub(super) fn reserve_bytes(capacity: u32, max_text_bytes_per_row: usize) -> usize {
    (capacity as usize)
        .saturating_mul(RETAINED_BYTES_PER_ROW.saturating_add(max_text_bytes_per_row))
        .saturating_add(DICTIONARY_HEADROOM)
}

pub(super) fn range_bytes(range: &std::ops::Range<u64>) -> Result<usize, StructureBatchError> {
    usize::try_from(range.end.saturating_sub(range.start)).map_err(|_| address_overflow())
}

pub(super) fn address_overflow() -> StructureBatchError {
    Diagnostic::new(Code::E1903)
        .with_message("BinaryCIF batch address exceeds the local address space")
        .into()
}

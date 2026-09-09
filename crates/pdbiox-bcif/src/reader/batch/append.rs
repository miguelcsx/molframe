//! Row filtering and direct append into the columnar structure builder.

use super::decoded::DecodedChunk;
use super::{ColumnStream, FIELD_COUNT, RowView};
use pdbiox_core::{ReadOptions, StructureBatchBuilder, StructureBatchError};

pub(super) fn append_rows(
    builder: &mut StructureBatchBuilder,
    columns: &[ColumnStream],
    fields: &[Option<usize>; FIELD_COUNT],
    chunks: &[DecodedChunk],
    options: &ReadOptions,
    rows: usize,
) -> Result<u32, StructureBatchError> {
    let mut emitted = 0_u32;
    for row in 0..rows {
        let view = RowView {
            columns,
            fields,
            chunks,
            row,
        };
        if view.accepted(options) {
            builder.push(view.record())?;
            emitted = emitted.saturating_add(1);
        }
    }
    Ok(emitted)
}

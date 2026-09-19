//! Arrow buffers retaining the immutable structure snapshot they alias.

use arrow::buffer::{Buffer, ScalarBuffer};
use arrow::datatypes::ArrowNativeType;
use arrow::error::{ArrowError, Result};
use molframe_core::{Structure, SymbolId};
use std::ptr::NonNull;
use std::sync::Arc;

#[derive(Debug)]
struct SnapshotAllocation {
    _structure: Structure,
}

/// Shared owner for every zero-copy buffer exported from one table snapshot.
#[derive(Clone, Debug)]
pub(crate) struct SnapshotOwner(Arc<SnapshotAllocation>);

impl SnapshotOwner {
    pub(crate) fn new(structure: &Structure) -> Self {
        Self(Arc::new(SnapshotAllocation {
            _structure: structure.clone(),
        }))
    }

    fn allocation(&self) -> Arc<dyn arrow::alloc::Allocation> {
        self.0.clone()
    }
}

pub(crate) fn f32_buffer(values: &[f32], owner: &SnapshotOwner) -> Result<ScalarBuffer<f32>> {
    typed_buffer(values, owner)
}

pub(crate) fn symbol_buffer(
    values: &[SymbolId],
    owner: &SnapshotOwner,
) -> Result<ScalarBuffer<u32>> {
    let bytes = values
        .len()
        .checked_mul(std::mem::size_of::<u32>())
        .ok_or_else(|| ArrowError::MemoryError("symbol buffer length overflow".to_owned()))?;
    let buffer = custom_buffer(values.as_ptr().cast::<u8>(), bytes, owner)?;
    Ok(ScalarBuffer::new(buffer, 0, values.len()))
}

fn typed_buffer<T: ArrowNativeType>(
    values: &[T],
    owner: &SnapshotOwner,
) -> Result<ScalarBuffer<T>> {
    let bytes = values
        .len()
        .checked_mul(std::mem::size_of::<T>())
        .ok_or_else(|| ArrowError::MemoryError("Arrow buffer length overflow".to_owned()))?;
    let buffer = custom_buffer(values.as_ptr().cast::<u8>(), bytes, owner)?;
    Ok(ScalarBuffer::new(buffer, 0, values.len()))
}

fn custom_buffer(pointer: *const u8, bytes: usize, owner: &SnapshotOwner) -> Result<Buffer> {
    let pointer = NonNull::new(pointer.cast_mut())
        .ok_or_else(|| ArrowError::MemoryError("null molframe buffer".to_owned()))?;
    let owner = owner.allocation();
    // SAFETY: the pointer and byte length come from a live immutable slice.
    // The shared snapshot owner retains the Structure allocation until Arrow
    // drops the final buffer, and public snapshots never mutate in place.
    Ok(unsafe { Buffer::from_custom_allocation(pointer, bytes, owner) })
}

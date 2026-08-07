//! Arrow buffers retaining the immutable structure snapshot they alias.

use arrow::buffer::{Buffer, ScalarBuffer};
use arrow::datatypes::ArrowNativeType;
use arrow::error::{ArrowError, Result};
use pdbiox_core::{Structure, SymbolId};
use std::ptr::NonNull;
use std::sync::Arc;

#[derive(Debug)]
struct SnapshotOwner {
    _structure: Structure,
}

pub(crate) fn f32_buffer(values: &[f32], owner: &Structure) -> Result<ScalarBuffer<f32>> {
    typed_buffer(values, owner)
}

pub(crate) fn symbol_buffer(values: &[SymbolId], owner: &Structure) -> Result<ScalarBuffer<u32>> {
    let bytes = values
        .len()
        .checked_mul(std::mem::size_of::<u32>())
        .ok_or_else(|| ArrowError::MemoryError("symbol buffer length overflow".to_owned()))?;
    let buffer = custom_buffer(values.as_ptr().cast::<u8>(), bytes, owner)?;
    Ok(ScalarBuffer::new(buffer, 0, values.len()))
}

fn typed_buffer<T: ArrowNativeType>(values: &[T], owner: &Structure) -> Result<ScalarBuffer<T>> {
    let bytes = values
        .len()
        .checked_mul(std::mem::size_of::<T>())
        .ok_or_else(|| ArrowError::MemoryError("Arrow buffer length overflow".to_owned()))?;
    let buffer = custom_buffer(values.as_ptr().cast::<u8>(), bytes, owner)?;
    Ok(ScalarBuffer::new(buffer, 0, values.len()))
}

fn custom_buffer(pointer: *const u8, bytes: usize, owner: &Structure) -> Result<Buffer> {
    let pointer = NonNull::new(pointer.cast_mut())
        .ok_or_else(|| ArrowError::MemoryError("null pdbiox buffer".to_owned()))?;
    let owner: Arc<dyn arrow::alloc::Allocation> = Arc::new(SnapshotOwner {
        _structure: owner.clone(),
    });
    // SAFETY: the pointer and byte length come from a live immutable slice.
    // `SnapshotOwner` retains the Structure allocation until Arrow drops the
    // final buffer, and public structure snapshots never mutate in place.
    Ok(unsafe { Buffer::from_custom_allocation(pointer, bytes, owner) })
}

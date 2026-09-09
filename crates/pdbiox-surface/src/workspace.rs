//! Complete byte accounting and fallible allocation for voxel surfaces.

use core::mem::size_of;

use crate::{Cavity, SasaError, SurfaceTriangle};

/// Checks the input plus promoted radii before the first operation-owned
/// allocation is requested.
pub(crate) fn ensure_surface_input(atom_count: usize, limit: usize) -> Result<(), SasaError> {
    ensure(
        limit,
        sum(&[input_bytes(atom_count)?, bytes_for::<f64>(atom_count)?])?,
    )
}

/// Checks all non-output phases of a cavity calculation before allocating its
/// cell state or flood stack.
pub(crate) fn ensure_cavity_grid(
    atom_count: usize,
    cell_count: usize,
    limit: usize,
) -> Result<(), SasaError> {
    let input = input_bytes(atom_count)?;
    let expanded = bytes_for::<f64>(atom_count)?;
    let state = bytes_for::<u8>(cell_count)?;
    let stack = bytes_for::<u32>(cell_count)?;
    let initial = sum(&[input, expanded, state])?;

    // No more than every other grid cell can be a separate six-connected
    // component. Doubling the output allowance also covers Vec reallocation.
    let component_count = cell_count
        .checked_add(1)
        .ok_or(SasaError::GridDimensionsOverflow)?
        / 2;
    let output = bytes_for::<Cavity>(component_count)?
        .checked_mul(2)
        .ok_or(SasaError::GridDimensionsOverflow)?;
    let flooded = sum(&[input, state, stack, output])?;

    ensure(limit, initial.max(flooded))
}

/// Checks all grid-only phases of an SES calculation before allocating its
/// cell state or flood stack.
pub(crate) fn ensure_ses_grid(
    atom_count: usize,
    cell_count: usize,
    longest_line: usize,
    limit: usize,
) -> Result<(), SasaError> {
    let input = input_bytes(atom_count)?;
    let expanded = bytes_for::<f64>(atom_count)?;
    let state = bytes_for::<u8>(cell_count)?;
    let stack = bytes_for::<u32>(cell_count)?;
    let distance = bytes_for::<f32>(cell_count)?;
    let distance_scratch = sum(&[
        bytes_for::<u32>(longest_line)?,
        bytes_for::<f32>(longest_line)?,
        bytes_for::<f64>(longest_line)?,
    ])?;
    let initial = sum(&[input, expanded, state])?;
    let flooded = sum(&[input, state, stack])?;
    let transformed = sum(&[input, state, distance])?;
    let distance_transform = sum(&[input, distance, distance_scratch])?;

    ensure(
        limit,
        initial
            .max(flooded)
            .max(transformed)
            .max(distance_transform),
    )
}

/// Checks the exact reserved SES output capacity together with the live
/// distance field before the triangle allocation occurs.
pub(crate) fn ensure_ses_output(
    atom_count: usize,
    cell_count: usize,
    triangle_capacity: usize,
    limit: usize,
) -> Result<(), SasaError> {
    let bytes = sum(&[
        input_bytes(atom_count)?,
        bytes_for::<f32>(cell_count)?,
        bytes_for::<SurfaceTriangle>(triangle_capacity)?,
    ])?;
    ensure(limit, bytes)
}

/// Allocates capacity without invoking the infallible global-OOM path.
pub(crate) fn empty_with_capacity<T>(capacity: usize) -> Result<Vec<T>, SasaError> {
    let bytes = bytes_for::<T>(capacity)?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(capacity)
        .map_err(|_| SasaError::AllocationFailed { bytes })?;
    Ok(values)
}

/// Allocates and initializes a bounded cell buffer fallibly.
pub(crate) fn filled<T: Clone>(length: usize, value: T) -> Result<Vec<T>, SasaError> {
    let mut values = empty_with_capacity(length)?;
    values.resize(length, value);
    Ok(values)
}

/// Returns the byte count of `count` values or a deterministic size error.
fn bytes_for<T>(count: usize) -> Result<usize, SasaError> {
    count
        .checked_mul(size_of::<T>())
        .ok_or(SasaError::GridDimensionsOverflow)
}

/// Includes caller-owned coordinates and radii because they remain resident
/// throughout every surface operation.
fn input_bytes(atom_count: usize) -> Result<usize, SasaError> {
    sum(&[
        bytes_for::<[f32; 3]>(atom_count)?,
        bytes_for::<f32>(atom_count)?,
    ])
}

fn sum(parts: &[usize]) -> Result<usize, SasaError> {
    let mut total = 0usize;
    for &part in parts {
        total = total
            .checked_add(part)
            .ok_or(SasaError::GridDimensionsOverflow)?;
    }
    Ok(total)
}

fn ensure(limit: usize, bytes: usize) -> Result<(), SasaError> {
    if bytes > limit {
        Err(SasaError::WorkspaceTooLarge { bytes, limit })
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "workspace_tests.rs"]
mod tests;

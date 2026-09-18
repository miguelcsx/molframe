//! Direct diagonalisation of the anisotropic Hessian.
//!
//! The Hessian is recovered one column at a time through the same matrix-free
//! kernel the operator uses, then diagonalised outright. That costs
//! `O((3N)^2)` bytes and `O((3N)^3)` time, which puts a single protein of a few
//! hundred residues in the low seconds and a few tens of megabytes, and refuses
//! anything larger against the caller's own ceiling rather than trying and
//! failing slowly.
//!
//! A restarted Krylov method would lift that ceiling in principle, and was
//! tried: an elastic network's slow modes are so tightly clustered that the
//! available implementation's re-orthogonalisation stops making progress and
//! spins. Lifting the ceiling properly needs a spectral transformation — a
//! sparse Cholesky of the shifted Hessian, which is nearly banded when sites
//! are ordered along the chain — rather than a better-tuned iteration. Until
//! that exists, refusing loudly beats hanging.

use nalgebra::{DMatrix, SymmetricEigen};

use super::hessian::Hessian;
use super::{AnisotropicNetworkModel, AnmError, AnmOptions, budget, canonical_mode};
use crate::network;

/// Working copies of the matrix the dense solve holds at once: the operator
/// itself and the eigenvector basis it is decomposed into.
const DENSE_MATRIX_COUNT: usize = 2;

pub(super) fn solve(
    hessian: &Hessian<'_>,
    sites: Vec<u32>,
    options: AnmOptions,
) -> Result<AnisotropicNetworkModel, AnmError> {
    let dimension = hessian.dimension();
    let matrix_bytes = checked_product(
        &[dimension, dimension, DENSE_MATRIX_COUNT, size_of::<f64>()],
        options,
    )?;
    let output_bytes =
        checked_product(&[dimension, options.mode_count, size_of::<f64>()], options)?;
    let site_bytes = checked_product(&[sites.capacity(), size_of::<u32>()], options)?;
    check_memory(
        checked_sum(
            &[
                hessian.owned_bytes(),
                site_bytes,
                matrix_bytes,
                output_bytes,
            ],
            options,
        )?,
        options,
    )?;

    // One application per basis vector recovers the operator column by column,
    // which keeps the only description of the Hessian the edge kernel rather
    // than a second transcription that could drift from it.
    let mut dense = DMatrix::zeros(dimension, dimension);
    let mut unit = vec![0.0_f64; dimension];
    let mut column = vec![0.0_f64; dimension];
    for index in 0..dimension {
        unit.fill(0.0);
        unit[index] = 1.0;
        column.fill(0.0);
        hessian.accumulate(&mut column, &unit, 1.0);
        for (row, value) in column.iter().enumerate() {
            dense[(row, index)] = *value;
        }
    }

    let eigen = SymmetricEigen::new(dense);
    let mut indexed: Vec<(f64, usize)> = eigen
        .eigenvalues
        .iter()
        .copied()
        .enumerate()
        .map(|(index, value)| (value, index))
        .collect();
    // A stable key: equal eigenvalues keep their column order, so the same
    // input reports the same modes in the same order on every run.
    indexed.sort_by(|left, right| left.0.total_cmp(&right.0).then(left.1.cmp(&right.1)));

    let zero_modes = indexed
        .iter()
        .take_while(|(value, _)| value.abs() <= options.zero_mode_tolerance)
        .count();
    let available = indexed.len() - zero_modes;
    if available < options.mode_count {
        return Err(AnmError::InsufficientModes {
            requested: options.mode_count,
            available,
        });
    }

    let chosen = &indexed[zero_modes..zero_modes + options.mode_count];
    Ok(AnisotropicNetworkModel {
        sites,
        eigenvalues: chosen.iter().map(|(value, _)| *value).collect(),
        modes: chosen
            .iter()
            .map(|(_, column)| {
                canonical_mode(displacements(
                    eigen.eigenvectors.column(*column).iter().copied(),
                ))
            })
            .collect(),
        zero_modes,
    })
}

/// Splits a flat `3N` vector into one displacement per site.
fn displacements(values: impl Iterator<Item = f64>) -> Vec<[f64; 3]> {
    let mut displacements = Vec::new();
    let mut current = [0.0_f64; 3];
    for (index, value) in values.enumerate() {
        current[index % 3] = value;
        if index % 3 == 2 {
            displacements.push(current);
        }
    }
    displacements
}

fn checked_product(factors: &[usize], options: AnmOptions) -> Result<usize, AnmError> {
    network::checked_product(factors, budget(options)).map_err(AnmError::from)
}

fn checked_sum(values: &[usize], options: AnmOptions) -> Result<usize, AnmError> {
    network::checked_sum(values, budget(options)).map_err(AnmError::from)
}

fn check_memory(required: usize, options: AnmOptions) -> Result<(), AnmError> {
    network::check_memory(required, budget(options)).map_err(AnmError::from)
}

//! Gaussian network model over explicitly selected interaction sites.

use nalgebra::{DMatrix, SymmetricEigen};
use pdbiox_core::selection::AtomSelection;
use pdbiox_spatial::{PeriodicBox, SpatialBackend, SpatialError, pairs_within};
use std::collections::BTreeMap;

/// Explicit Gaussian-network construction and solve policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GnmOptions {
    /// Maximum separation defining a network spring.
    pub contact_distance: f32,
    /// Number of non-zero modes requested.
    pub mode_count: usize,
    /// Eigenvalues at or below this magnitude are classified as zero modes.
    pub zero_mode_tolerance: f64,
    /// Maximum bytes allowed for the dense Kirchhoff matrix.
    pub memory_limit_bytes: usize,
    /// Spatial implementation used to build contacts.
    pub backend: SpatialBackend,
}

/// Non-zero Gaussian network modes in increasing eigenvalue order.
#[derive(Clone, Debug, PartialEq)]
pub struct GaussianNetworkModel {
    /// Selected site indices defining vector element order.
    pub sites: Vec<u32>,
    /// Non-zero eigenvalues in increasing order.
    pub eigenvalues: Vec<f64>,
    /// Corresponding canonical-sign eigenvectors.
    pub modes: Vec<Vec<f64>>,
    /// Number of zero modes, including disconnected-component translations.
    pub zero_modes: usize,
}

/// Why a Gaussian network could not be constructed or solved.
#[derive(Debug, thiserror::Error)]
pub enum GnmError {
    /// Network controls are invalid.
    #[error("GNM cutoff, mode count, tolerance, and memory limit must be valid and explicit")]
    InvalidOptions,
    /// At least two selected sites are required.
    #[error("GNM requires at least two selected sites")]
    TooFewSites,
    /// Dense matrix allocation exceeds the caller's ceiling.
    #[error("GNM matrix requires {required} bytes, over the {limit} byte limit")]
    MemoryLimit {
        /// Required matrix bytes.
        required: usize,
        /// Caller-provided ceiling.
        limit: usize,
    },
    /// Fewer non-zero modes exist than requested.
    #[error("requested {requested} non-zero modes, but only {available} exist")]
    InsufficientModes {
        /// Requested non-zero modes.
        requested: usize,
        /// Available non-zero modes.
        available: usize,
    },
    /// Spatial network construction failed.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
}

/// Constructs and diagonalizes a Gaussian network Kirchhoff matrix.
///
/// Every contact receives the conventional unit spring weight. Chemical site
/// choice, contact distance, zero-mode tolerance, memory ceiling, and periodic
/// geometry are all explicit inputs.
///
/// # Errors
///
/// Returns [`GnmError`] for invalid controls, selections, memory, or mode count.
pub fn gaussian_network_model(
    positions: &[[f32; 3]],
    sites: &AtomSelection,
    options: GnmOptions,
    periodic: Option<&PeriodicBox>,
) -> Result<GaussianNetworkModel, GnmError> {
    validate_options(options)?;
    let selected: Vec<u32> = sites.into_iter().collect();
    if selected.len() < 2 {
        return Err(GnmError::TooFewSites);
    }
    let required = selected
        .len()
        .checked_mul(selected.len())
        .and_then(|elements| elements.checked_mul(size_of::<f64>()))
        .ok_or(GnmError::MemoryLimit {
            required: usize::MAX,
            limit: options.memory_limit_bytes,
        })?;
    if required > options.memory_limit_bytes {
        return Err(GnmError::MemoryLimit {
            required,
            limit: options.memory_limit_bytes,
        });
    }
    let contacts = pairs_within(
        positions,
        sites,
        sites,
        options.contact_distance,
        options.backend,
        periodic,
    )?;
    let lookup: BTreeMap<u32, usize> = selected
        .iter()
        .copied()
        .enumerate()
        .map(|(local, atom)| (atom, local))
        .collect();
    let mut kirchhoff = DMatrix::zeros(selected.len(), selected.len());
    for pair in contacts {
        let Some(&left) = lookup.get(&pair.first) else {
            continue;
        };
        let Some(&right) = lookup.get(&pair.second) else {
            continue;
        };
        kirchhoff[(left, right)] = -1.0;
        kirchhoff[(right, left)] = -1.0;
        kirchhoff[(left, left)] += 1.0;
        kirchhoff[(right, right)] += 1.0;
    }
    modes_from_matrix(kirchhoff, selected, options)
}

fn validate_options(options: GnmOptions) -> Result<(), GnmError> {
    if options.contact_distance.is_finite()
        && options.contact_distance > 0.0
        && options.mode_count > 0
        && options.zero_mode_tolerance.is_finite()
        && options.zero_mode_tolerance >= 0.0
        && options.memory_limit_bytes > 0
    {
        Ok(())
    } else {
        Err(GnmError::InvalidOptions)
    }
}

fn modes_from_matrix(
    kirchhoff: DMatrix<f64>,
    sites: Vec<u32>,
    options: GnmOptions,
) -> Result<GaussianNetworkModel, GnmError> {
    let eigen = SymmetricEigen::new(kirchhoff);
    let mut indexed: Vec<(f64, usize)> = eigen
        .eigenvalues
        .iter()
        .copied()
        .enumerate()
        .map(|(index, value)| (value, index))
        .collect();
    indexed.sort_by(|left, right| left.0.total_cmp(&right.0));
    let zero_modes = indexed
        .iter()
        .take_while(|(value, _)| value.abs() <= options.zero_mode_tolerance)
        .count();
    let available = indexed.len() - zero_modes;
    if available < options.mode_count {
        return Err(GnmError::InsufficientModes {
            requested: options.mode_count,
            available,
        });
    }
    let chosen = &indexed[zero_modes..zero_modes + options.mode_count];
    let eigenvalues = chosen.iter().map(|(value, _)| *value).collect();
    let modes = chosen
        .iter()
        .map(|(_, column)| {
            canonical_mode(eigen.eigenvectors.column(*column).iter().copied().collect())
        })
        .collect();
    Ok(GaussianNetworkModel {
        sites,
        eigenvalues,
        modes,
        zero_modes,
    })
}

fn canonical_mode(mut mode: Vec<f64>) -> Vec<f64> {
    if mode
        .iter()
        .copied()
        .find(|value| *value != 0.0)
        .is_some_and(|value| value < 0.0)
    {
        for value in &mut mode {
            *value = -*value;
        }
    }
    mode
}

#[cfg(test)]
#[path = "gnm_tests.rs"]
mod tests;

//! Gaussian network model over explicitly selected interaction sites.
//!
//! Contacts and every matrix-vector product cost `O(N + E)`, where `N` is the
//! number of selected sites and `E` is the number of cutoff contacts. Only the
//! small-site micro-kernel materialises an `N x N` matrix.

#[path = "gnm/solver.rs"]
mod solver;

use crate::network::{ContactGraph, NetworkBudget, NetworkError};
use pdbiox_core::ExecutionContext;
use pdbiox_core::parallel::ReductionPolicy;
use pdbiox_core::selection::AtomSelection;
use pdbiox_spatial::{PeriodicBox, SpatialBackend, SpatialError};

/// Explicit Gaussian-network construction and solve policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GnmOptions {
    /// Maximum separation defining a network spring.
    pub contact_distance: f32,
    /// Number of non-zero modes requested.
    pub mode_count: usize,
    /// Eigenvalues at or below this magnitude are classified as zero modes.
    pub zero_mode_tolerance: f64,
    /// Maximum bytes for graph construction, eigensolver workspace, and output.
    pub memory_limit_bytes: usize,
    /// Spatial implementation used to build contacts.
    pub backend: SpatialBackend,
    /// Whether the eigensolver may reorder its floating-point reductions.
    ///
    /// `Deterministic` — the default — solves sequentially, so the eigenvalues
    /// are bit-identical on every run and machine (FR-515). `Fast` lets the
    /// solver use the execution context's native pool and returns results that
    /// agree to within the solver tolerance rather than exactly.
    pub reduction: ReductionPolicy,
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
    /// Sparse graph, solver workspace, and output exceed the caller's ceiling.
    #[error("GNM workspace requires {required} bytes, over the {limit} byte limit")]
    MemoryLimit {
        /// Required workspace bytes.
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
    /// The restarted sparse solve exhausted its convergence budget.
    #[error("GNM sparse eigensolver converged {converged} of {requested} requested modes")]
    Convergence {
        /// Modes requested from the sparse solver.
        requested: usize,
        /// Modes that reached the residual threshold.
        converged: usize,
    },
    /// Spatial network construction failed.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
}

/// Constructs and diagonalizes a Gaussian network Kirchhoff operator.
///
/// Every contact receives the conventional unit spring weight. Chemical site
/// choice, contact distance, zero-mode tolerance, memory ceiling, and periodic
/// geometry are all explicit inputs.
///
/// # Errors
///
/// Returns [`GnmError`] for invalid controls, selections, memory, convergence,
/// or mode count.
pub fn gaussian_network_model(
    positions: &[[f32; 3]],
    sites: &AtomSelection,
    options: GnmOptions,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<GaussianNetworkModel, GnmError> {
    validate_options(options)?;
    let site_count = usize::try_from(sites.len()).map_err(|_| GnmError::MemoryLimit {
        required: usize::MAX,
        limit: options.memory_limit_bytes,
    })?;
    if site_count < 2 {
        return Err(GnmError::TooFewSites);
    }
    let site_bytes = solver::checked_product(&[site_count, size_of::<u32>()], options)?;
    solver::check_memory(site_bytes, options)?;
    let selected: Vec<u32> = sites.into_iter().collect();
    validate_site_indices(&selected, positions.len())?;
    let graph = ContactGraph::build(
        positions,
        sites,
        &selected,
        budget(options),
        periodic,
        context,
    )?;
    solver::solve(&graph, selected, options)
}

fn validate_site_indices(selected: &[u32], position_count: usize) -> Result<(), GnmError> {
    for &atom in selected {
        let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        if index >= position_count {
            return Err(SpatialError::AtomOutOfBounds(atom).into());
        }
    }
    Ok(())
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

pub(super) fn canonical_mode(mut mode: Vec<f64>) -> Vec<f64> {
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

/// The subset of the public options the shared contact build actually reads.
fn budget(options: GnmOptions) -> NetworkBudget {
    NetworkBudget {
        contact_distance: options.contact_distance,
        memory_limit_bytes: options.memory_limit_bytes,
        backend: options.backend,
    }
}

impl From<NetworkError> for GnmError {
    fn from(error: NetworkError) -> Self {
        match error {
            NetworkError::MemoryLimit { required, limit } => Self::MemoryLimit { required, limit },
            NetworkError::Spatial(error) => Self::Spatial(error),
        }
    }
}

#[cfg(test)]
#[path = "gnm_tests.rs"]
mod tests;

//! Anisotropic network model over explicitly selected interaction sites.
//!
//! The Gaussian model answers *how much* each site moves; this one answers
//! *which way*. Each cutoff contact contributes a rank-one `3x3` block along
//! the bond direction, so a mode is a field of three-component displacements
//! rather than a scalar per site — which is what a porcupine plot, a normal
//! mode animation, or a physically grounded apo-to-holo morph actually needs.
//!
//! Contacts cost `O(N + E)` to build and one application of the Hessian costs
//! `O(N + E)`, with the bond directions normalised once rather than per apply.
//! The eigenproblem itself is solved directly, which costs `O((3N)^2)` bytes
//! and `O((3N)^3)` time. Measured on one core: a hundred sites in about thirty
//! milliseconds, two hundred and fifty in half a second, five hundred in three
//! and a half. That covers a single domain or a small complex comfortably;
//! beyond it the caller's memory ceiling refuses the solve and names the bytes
//! it would have needed, rather than starting one that will not finish.

#[path = "anm/hessian.rs"]
mod hessian;
#[path = "anm/solver.rs"]
mod solver;

use crate::network::{ContactGraph, NetworkBudget, NetworkError};
use crate::numeric::f64_to_f32;
use hessian::Hessian;
use molframe_core::ExecutionContext;
use molframe_core::selection::AtomSelection;
use molframe_spatial::{PeriodicBox, SpatialBackend, SpatialError};

/// Explicit anisotropic-network construction and solve policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnmOptions {
    /// Maximum separation defining a network spring.
    ///
    /// Alpha-carbon networks are conventionally sprung well beyond a covalent
    /// distance — the spring stands in for every interaction that keeps two
    /// residues near each other, not for a bond.
    pub contact_distance: f32,
    /// Number of non-zero modes requested.
    pub mode_count: usize,
    /// Eigenvalues at or below this magnitude are classified as zero modes.
    pub zero_mode_tolerance: f64,
    /// Maximum bytes for graph construction, eigensolver workspace, and output.
    pub memory_limit_bytes: usize,
    /// Spatial implementation used to build contacts.
    pub backend: SpatialBackend,
}

/// Non-zero anisotropic network modes in increasing eigenvalue order.
#[derive(Clone, Debug, PartialEq)]
pub struct AnisotropicNetworkModel {
    /// Selected site indices defining displacement element order.
    pub sites: Vec<u32>,
    /// Non-zero eigenvalues in increasing order.
    ///
    /// A mode's eigenvalue is its stiffness: the lowest ones are the slow,
    /// collective motions a structure actually performs.
    pub eigenvalues: Vec<f64>,
    /// Per mode, one unit-normalised displacement per site.
    pub modes: Vec<Vec<[f64; 3]>>,
    /// Directions with no restoring force.
    ///
    /// Six for a connected structure with extent about every axis — three
    /// translations and three rotations. More means the network is degenerate:
    /// a collinear chain, for instance, has springs that resist nothing across
    /// their own axis, and this count says so rather than hiding it.
    pub zero_modes: usize,
}

impl AnisotropicNetworkModel {
    /// Mean-square fluctuation of every site, in the model's arbitrary units.
    ///
    /// Each mode contributes its squared displacement scaled by the inverse of
    /// its stiffness, so the slow modes dominate exactly as they do in a real
    /// temperature factor. `O(modes x sites)`.
    #[must_use]
    pub fn fluctuations(&self) -> Vec<f64> {
        let mut fluctuations = vec![0.0; self.sites.len()];
        for (eigenvalue, mode) in self.eigenvalues.iter().zip(&self.modes) {
            if *eigenvalue <= 0.0 {
                continue;
            }
            for (fluctuation, displacement) in fluctuations.iter_mut().zip(mode) {
                let squared = displacement.iter().map(|value| value * value).sum::<f64>();
                *fluctuation += squared / eigenvalue;
            }
        }
        fluctuations
    }

    /// Overlap of each mode with a target displacement field.
    ///
    /// Feeding the apo-to-holo difference in gives the amplitude to travel
    /// along every mode; scaling those amplitudes by an animation parameter and
    /// summing them traces a path that respects the network's own stiffness
    /// instead of dragging atoms through each other in straight lines.
    /// `O(modes x sites)`.
    ///
    /// The displacement field is indexed like [`Self::sites`]; a shorter slice
    /// yields no overlaps.
    #[must_use]
    pub fn project(&self, displacement: &[[f64; 3]]) -> Vec<f64> {
        if displacement.len() < self.sites.len() {
            return Vec::new();
        }
        self.modes
            .iter()
            .map(|mode| {
                mode.iter()
                    .zip(displacement)
                    .map(|(basis, target)| {
                        basis[0] * target[0] + basis[1] * target[1] + basis[2] * target[2]
                    })
                    .sum()
            })
            .collect()
    }

    /// Writes `positions + sum(amplitude x mode)` into `out`.
    ///
    /// This is the animation step: one call per frame with amplitudes scaled by
    /// the frame's phase produces a deformation that never leaves the subspace
    /// the network permits. Amplitudes past the available modes are ignored, as
    /// are sites past the shorter of `positions` and `out`. `O(modes x sites)`.
    pub fn displace(&self, positions: &[[f32; 3]], amplitudes: &[f64], out: &mut [[f32; 3]]) {
        let count = self.sites.len().min(positions.len()).min(out.len());
        out[..count].copy_from_slice(&positions[..count]);
        for (amplitude, mode) in amplitudes.iter().zip(&self.modes) {
            if *amplitude == 0.0 {
                continue;
            }
            for (position, displacement) in out[..count].iter_mut().zip(mode) {
                for (axis, offset) in position.iter_mut().zip(displacement) {
                    // The accumulation is in f64 and only the result narrows,
                    // so a long sum of small amplitudes does not drift.
                    *axis = f64_to_f32(f64::from(*axis) + amplitude * offset);
                }
            }
        }
    }
}

/// Why an anisotropic network could not be constructed or solved.
#[derive(Debug, thiserror::Error)]
pub enum AnmError {
    /// Network controls are invalid.
    #[error("ANM cutoff, mode count, tolerance, and memory limit must be valid and explicit")]
    InvalidOptions,
    /// At least two selected sites are required.
    #[error("ANM requires at least two selected sites")]
    TooFewSites,
    /// Sparse graph, solver workspace, and output exceed the caller's ceiling.
    #[error("ANM workspace requires {required} bytes, over the {limit} byte limit")]
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
    #[error("ANM sparse eigensolver converged {converged} of {requested} requested modes")]
    Convergence {
        /// Modes requested from the sparse solver.
        requested: usize,
        /// Modes that reached the residual threshold.
        converged: usize,
    },
    /// Two selected sites occupy the same point, leaving no bond direction.
    #[error("ANM sites {first} and {second} are contacts at zero separation")]
    CoincidentSites {
        /// First site index.
        first: u32,
        /// Second site index.
        second: u32,
    },
    /// Spatial network construction failed.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
}

/// Constructs and diagonalizes an anisotropic network Hessian.
///
/// Every contact receives the conventional unit spring weight. Chemical site
/// choice, contact distance, zero-mode tolerance, memory ceiling, and periodic
/// geometry are all explicit inputs.
///
/// # Errors
///
/// Returns [`AnmError`] for invalid controls, selections, memory, convergence,
/// coincident sites, or mode count.
pub fn anisotropic_network_model(
    positions: &[[f32; 3]],
    sites: &AtomSelection,
    options: AnmOptions,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<AnisotropicNetworkModel, AnmError> {
    validate_options(options)?;
    let site_count = usize::try_from(sites.len()).map_err(|_| AnmError::MemoryLimit {
        required: usize::MAX,
        limit: options.memory_limit_bytes,
    })?;
    if site_count < 2 {
        return Err(AnmError::TooFewSites);
    }
    let budget = budget(options);
    let site_bytes = crate::network::checked_product(&[site_count, size_of::<u32>()], budget)?;
    crate::network::check_memory(site_bytes, budget)?;
    let selected: Vec<u32> = sites.into_iter().collect();
    validate_site_indices(&selected, positions.len())?;
    let graph = ContactGraph::build(positions, sites, &selected, budget, periodic, context)?;
    let hessian = Hessian::build(&graph, positions, &selected, budget)?;
    solver::solve(&hessian, selected, options)
}

fn validate_site_indices(selected: &[u32], position_count: usize) -> Result<(), AnmError> {
    for &atom in selected {
        let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        if index >= position_count {
            return Err(SpatialError::AtomOutOfBounds(atom).into());
        }
    }
    Ok(())
}

fn validate_options(options: AnmOptions) -> Result<(), AnmError> {
    if options.contact_distance.is_finite()
        && options.contact_distance > 0.0
        && options.mode_count > 0
        && options.zero_mode_tolerance.is_finite()
        && options.zero_mode_tolerance >= 0.0
        && options.memory_limit_bytes > 0
    {
        Ok(())
    } else {
        Err(AnmError::InvalidOptions)
    }
}

/// The subset of the public options the shared contact build actually reads.
pub(super) fn budget(options: AnmOptions) -> NetworkBudget {
    NetworkBudget {
        contact_distance: options.contact_distance,
        memory_limit_bytes: options.memory_limit_bytes,
        backend: options.backend,
    }
}

impl From<NetworkError> for AnmError {
    fn from(error: NetworkError) -> Self {
        match error {
            NetworkError::MemoryLimit { required, limit } => Self::MemoryLimit { required, limit },
            NetworkError::Spatial(error) => Self::Spatial(error),
        }
    }
}

/// Fixes the arbitrary sign an eigensolver returns, so the same input always
/// reports the same mode rather than one that flips between runs.
pub(super) fn canonical_mode(mut mode: Vec<[f64; 3]>) -> Vec<[f64; 3]> {
    let leading = mode.iter().flatten().copied().find(|value| *value != 0.0);
    if leading.is_some_and(|value| value < 0.0) {
        for displacement in &mut mode {
            for value in displacement {
                *value = -*value;
            }
        }
    }
    mode
}

#[cfg(test)]
#[path = "anm_tests.rs"]
mod tests;

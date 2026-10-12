//! Native molecular-surface providers with bounded files and cooperative shutdown.
//!
//! Input writing and output decoding are linear in atoms and output size. Solver
//! internals belong to the caller-selected executable; its private heap is not
//! charged to Rust's memory account. Time and temporary files are monitored every
//! 20 ms. The decoded mesh retains a conservative shared memory reservation.

mod coordinates;
mod output;
mod process;
mod provenance;

use crate::{MolecularSurfaceMesh, SurfaceOutputError, SurfaceOutputLimits};
use coordinates::NativeAtoms;
pub use coordinates::SurfaceCoordinateSample;
use molframe_core::execution::{ExecutionContext, MemoryReservation};
use std::{path::Path, time::Duration};

/// The actual scientific provider, rather than an approximation under its name.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NativeSurfaceProvider {
    /// Varshney's rolling-probe triangle generator.
    Surf,
    /// Sanner's analytical rolling-probe triangulation.
    Msms,
}

/// Borrowed coordinate and radius columns plus the selected original atom rows.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceAtoms<'a> {
    /// Full molecular coordinates; never copied into another molecular model.
    pub positions: &'a [[f32; 3]],
    /// One positive van der Waals radius for every coordinate row.
    pub radii: &'a [f32],
    /// Input sphere order, retaining each atom's original row identity.
    pub rows: &'a [u32],
}

/// Scientific controls and resource ceilings for one native provider invocation.
#[derive(Clone, Copy, Debug)]
pub struct NativeSurfaceOptions {
    /// Rolling solvent-probe radius in angstroms.
    pub probe_radius: f32,
    /// MSMS triangulation density in vertices per square angstrom.
    pub density: f32,
    /// Include disconnected MSMS surface components.
    pub all_components: bool,
    /// Maximum accepted MSMS components, including disconnected spheres and cavities.
    pub maximum_components: usize,
    /// Total private files, including input and diagnostic output.
    pub temporary_bytes: u64,
    /// Maximum elapsed solver runtime.
    pub timeout: Duration,
    /// Complete decoded output limits.
    pub output: SurfaceOutputLimits,
}

impl Default for NativeSurfaceOptions {
    fn default() -> Self {
        Self {
            probe_radius: 1.4,
            density: 3.0,
            all_components: true,
            maximum_components: 4096,
            temporary_bytes: 1_000_000_000,
            timeout: Duration::from_secs(120),
            output: SurfaceOutputLimits::default(),
        }
    }
}

/// Complete provider geometry retaining its execution-memory charge.
#[derive(Debug)]
pub struct NativeSurfaceMesh {
    /// Coordinates, analytic normals, triangles and original molecular atom rows.
    pub mesh: MolecularSurfaceMesh,
    _memory: MemoryReservation,
    provider: NativeSurfaceProvider,
    options: NativeSurfaceOptions,
}

impl NativeSurfaceMesh {
    /// Actual scientific provider that produced this complete geometry.
    #[must_use]
    pub const fn provider(&self) -> NativeSurfaceProvider {
        self.provider
    }
    /// Scientific controls and resource ceilings used for the computation.
    #[must_use]
    pub const fn options(&self) -> NativeSurfaceOptions {
        self.options
    }
}

/// Failure to obtain a complete surface from a native executable.
#[derive(Debug, thiserror::Error)]
pub enum NativeSurfaceError {
    /// Input columns or controls are inconsistent.
    #[error("invalid native surface input: {0}")]
    Invalid(&'static str),
    /// Native filesystem or process operation failed.
    #[error(transparent)]
    Io(#[from] std::io::Error),
    /// Shared temporary-storage ceiling was exceeded or storage is disabled.
    #[error(transparent)]
    Storage(#[from] molframe_core::execution::SpillError),
    /// Decoded geometry cannot fit the shared execution-memory ceiling.
    #[error(transparent)]
    Memory(#[from] molframe_core::execution::MemoryBudgetError),
    /// Cooperative cancellation was requested.
    #[error(transparent)]
    Cancelled(#[from] molframe_core::execution::Cancelled),
    /// The provider exceeded its elapsed-time ceiling.
    #[error("native surface provider exceeded its time limit")]
    Timeout,
    /// The provider returned unsuccessful process status and bounded diagnostic output.
    #[error("native surface provider failed with status {status}: {diagnostic}")]
    Failed {
        /// Native exit code or terminating signal.
        status: std::process::ExitStatus,
        /// Final diagnostic bytes decoded as text, limited to eight kibibytes.
        diagnostic: String,
    },
    /// The provider output was malformed, incomplete or over budget.
    #[error(transparent)]
    Output(#[from] SurfaceOutputError),
}

/// Runs an explicitly supplied executable; no downloading or executable discovery.
/// A private working directory is removed on every success or failure path.
///
/// # Errors
/// Refuses invalid atom columns, controls, unsupported components, exhausted
/// resources, cancellation, unsuccessful providers or invalid geometry.
pub fn native_molecular_surface(
    executable: &Path,
    provider: NativeSurfaceProvider,
    atoms: SurfaceAtoms<'_>,
    options: NativeSurfaceOptions,
    context: &ExecutionContext,
) -> Result<NativeSurfaceMesh, NativeSurfaceError> {
    run(
        executable,
        provider,
        NativeAtoms {
            atoms,
            sample: None,
        },
        options,
        context,
    )
}

/// Computes a native analytical surface for a borrowed linear trajectory sample.
/// Sampling occurs only for selected input spheres and original atom identities
/// survive output decoding. Neither coordinate frame is copied or materialized.
///
/// # Errors
/// Rejects unmatched frames, invalid weights or any ordinary native provider failure.
pub fn native_molecular_surface_at_sample(
    executable: &Path,
    provider: NativeSurfaceProvider,
    atoms: SurfaceAtoms<'_>,
    sample: SurfaceCoordinateSample<'_>,
    options: NativeSurfaceOptions,
    context: &ExecutionContext,
) -> Result<NativeSurfaceMesh, NativeSurfaceError> {
    run(
        executable,
        provider,
        NativeAtoms {
            atoms,
            sample: Some(sample),
        },
        options,
        context,
    )
}

fn run(
    executable: &Path,
    provider: NativeSurfaceProvider,
    atoms: NativeAtoms<'_>,
    options: NativeSurfaceOptions,
    context: &ExecutionContext,
) -> Result<NativeSurfaceMesh, NativeSurfaceError> {
    validate(atoms, options)?;
    context.cancellation().check()?;
    let executable = std::fs::canonicalize(executable)?;
    let maximum_files = options
        .maximum_components
        .checked_mul(2)
        .and_then(|value| value.checked_add(2))
        .ok_or(NativeSurfaceError::Invalid(
            "component limit exceeds file capacity",
        ))?;
    let workspace =
        context.temporary_workspace_with_file_limit(options.temporary_bytes, maximum_files)?;
    process::write_input(workspace.path(), provider, atoms, context)?;
    workspace.check_size()?;
    process::run(&executable, provider, options, &workspace, context)?;
    workspace.check_size()?;
    let (mesh, memory) = output::read(workspace.path(), provider, atoms, options, context)?;
    Ok(NativeSurfaceMesh {
        mesh,
        _memory: memory,
        provider,
        options,
    })
}

fn validate(
    atoms: NativeAtoms<'_>,
    options: NativeSurfaceOptions,
) -> Result<(), NativeSurfaceError> {
    atoms.validate()?;
    if atoms.atoms.rows.is_empty() || atoms.atoms.positions.len() != atoms.atoms.radii.len() {
        return Err(NativeSurfaceError::Invalid(
            "empty selection or unmatched columns",
        ));
    }
    if !options.probe_radius.is_finite()
        || options.probe_radius <= 0.0
        || !options.density.is_finite()
        || options.density <= 0.0
        || options.timeout.is_zero()
        || options.temporary_bytes == 0
        || options.maximum_components == 0
    {
        return Err(NativeSurfaceError::Invalid(
            "controls must be finite and positive",
        ));
    }
    for &row in atoms.atoms.rows {
        let Some(position) = atoms.position(row) else {
            return Err(NativeSurfaceError::Invalid("selected atom exceeds columns"));
        };
        let radius = atoms.atoms.radii[row as usize];
        if !position.iter().all(|value| value.is_finite()) || !radius.is_finite() || radius <= 0.0 {
            return Err(NativeSurfaceError::Invalid(
                "atom coordinates or radii are invalid",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "native_solver_tests.rs"]
mod tests;

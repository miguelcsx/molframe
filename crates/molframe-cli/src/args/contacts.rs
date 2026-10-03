//! Arguments for residue contact maps, native contacts and normal modes.

use clap::{Args, ValueEnum};
use std::path::PathBuf;

/// Residue pairs whose closest atoms lie within a cutoff.
#[derive(Args, Debug)]
pub(crate) struct ContactMapArguments {
    /// The file to read.
    pub(crate) input: PathBuf,
    /// Largest atom-to-atom distance counted as a contact, in Angstrom.
    #[arg(long)]
    pub(crate) cutoff: f32,
    /// Smallest residue-index separation of a kept pair; 0 keeps every pair.
    #[arg(long)]
    pub(crate) min_separation: u32,
}

/// The fraction of a reference's contacts that a target keeps.
#[derive(Args, Debug)]
pub(crate) struct NativeContactArguments {
    /// The structure whose contacts are native.
    pub(crate) reference: PathBuf,
    /// The structure that shares its atom numbering and is scored.
    pub(crate) target: PathBuf,
    /// Distance defining a reference contact, in Angstrom.
    #[arg(long)]
    pub(crate) cutoff: f32,
    /// A contact is kept when the target distance is at most this multiple of the cutoff (1.0 demands it be within the cutoff again).
    #[arg(long)]
    pub(crate) retention: f32,
}

/// The elastic network model.
#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum NetworkModel {
    /// Gaussian network model: isotropic fluctuations of one scalar per site.
    Gnm,
    /// Anisotropic network model: a displacement vector per site.
    Anm,
}

/// What the normal-mode analysis prints.
#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum ModeTable {
    /// One row per non-zero mode: its eigenvalue.
    Modes,
    /// One row per site: its mean-square fluctuation over the computed modes.
    Fluctuations,
}

/// Whether the eigensolver may reorder floating-point reductions.
#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum ReductionChoice {
    /// Bit-identical output at any worker count.
    Deterministic,
    /// Faster; results agree to within the solver tolerance.
    Fast,
}

impl From<ReductionChoice> for molframe_core::parallel::ReductionPolicy {
    fn from(value: ReductionChoice) -> Self {
        match value {
            ReductionChoice::Deterministic => Self::Deterministic,
            ReductionChoice::Fast => Self::Fast,
        }
    }
}

/// Elastic-network normal-mode analysis over selected sites.
#[derive(Args, Debug)]
pub(crate) struct NormalModeArguments {
    /// The file to read.
    pub(crate) input: PathBuf,
    /// Which elastic network to build.
    #[arg(long, value_enum)]
    pub(crate) network: NetworkModel,
    /// Selection of the network sites, such as the alpha carbons.
    #[arg(long)]
    pub(crate) sites: String,
    /// Largest separation that defines a network spring, in Angstrom.
    #[arg(long)]
    pub(crate) contact_distance: f32,
    /// Number of non-zero modes requested.
    #[arg(long)]
    pub(crate) mode_count: usize,
    /// Eigenvalues at or below this magnitude are zero modes.
    #[arg(long)]
    pub(crate) zero_mode_tolerance: f64,
    /// Byte ceiling of graph construction, solver workspace and output.
    #[arg(long)]
    pub(crate) memory_limit: usize,
    /// Reduction order of the GNM eigensolver (required for `gnm`).
    #[arg(long, value_enum)]
    pub(crate) reduction: Option<ReductionChoice>,
    /// What to print.
    #[arg(long, value_enum)]
    pub(crate) table: ModeTable,
}

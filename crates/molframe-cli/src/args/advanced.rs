//! Arguments for reproducibility, functional evaluation and batch workflows.

use super::workflows::CcdArguments;
use clap::{Args, Subcommand};
use std::path::PathBuf;

#[derive(Args, Debug)]
pub(crate) struct AuditArguments {
    pub(crate) input: PathBuf,
    /// Declarative selection whose membership is audited.
    pub(crate) query: String,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) assembly_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) model_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) altloc_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) identifier_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) missing_atom_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) hydrogen_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) atom_equivalence_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) symmetry_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) alignment_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) precision_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) periodic_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) vdw_radii_values: Vec<String>,
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) contact_values: Vec<String>,
    /// Tolerance alternatives as RELATIVE:ABSOLUTE.
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub(crate) float_tolerance_values: Vec<String>,
    /// Hard bound on Cartesian policy points.
    #[arg(long)]
    pub(crate) max_runs: usize,
}

#[derive(Subcommand, Debug)]
pub(crate) enum FxCommand {
    /// Evaluate a structure against a declarative motif and verdict profile.
    Evaluate {
        input: PathBuf,
        #[arg(long)]
        motif: PathBuf,
        #[command(flatten)]
        chemistry: CcdArguments,
        #[arg(long)]
        mapping_limit: usize,
        #[arg(long)]
        measurement_limit: usize,
        /// Relative convergence tolerance for planarity measurements.
        #[arg(long)]
        plane_relative_tolerance: f64,
        /// Hard ceiling on cyclic-Jacobi sweeps for planarity measurements.
        #[arg(long)]
        plane_maximum_sweeps: usize,
    },
}

#[derive(Subcommand, Debug)]
pub(crate) enum BatchCommand {
    /// Summarise every structure matching explicit paths or glob patterns.
    Info {
        #[arg(num_args = 1..)]
        inputs: Vec<String>,
    },
    /// Convert every matching structure to one explicitly named format.
    Convert {
        #[arg(num_args = 1..)]
        inputs: Vec<String>,
        #[arg(long)]
        to: String,
        #[arg(long)]
        outdir: PathBuf,
    },
}

//! The pdbiox command line.
//!
//! Argument parsing, output formatting and exit codes over the library. No
//! analysis happens here: anything worth doing is worth doing in the library,
//! where the Rust callers get it too.
//!
//! Results go to standard output and findings go to standard error, always, so
//! a pipeline is never contaminated by warnings.

#![forbid(unsafe_code)]

mod analysis_commands;
mod args;
mod audit_commands;
mod batch_commands;
mod chemistry;
mod choices;
mod commands;
mod comparison_commands;
mod diff_commands;
mod dispatch;
mod exit;
mod fx_commands;
mod global;
mod inspect;
mod intrinsic_commands;
mod man_commands;
mod network_commands;
mod report;
mod selection_commands;
mod sequence_commands;
mod system;
mod trajectory;
mod validation_commands;

use clap::{Parser, Subcommand};
use std::path::PathBuf;
use std::process::ExitCode;

pub(crate) use args::{
    AltlocArgument, AssemblyArgument, AuditArguments, BatchCommand, CcdArguments, EnsembleCommand,
    FxCommand, GeometryCommand, KmerOperation, MatrixChoice, ModelArgument, NamespaceArgument,
    PairwiseMode, SequenceCommand, SequenceFormat, SurfaceArguments, SystemCommand,
    TrajectoryCommand, TreeMethod,
};
pub(crate) use choices::{
    CompletionShell, EmptyLddtChoice, MetricChoice, RadiusChoice, ValidationChoice,
};

/// A batteries-included structural bioinformatics engine.
#[derive(Parser, Debug)]
#[command(name = "pdbiox", version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,

    #[command(flatten)]
    global: global::GlobalOptions,
}

#[derive(Subcommand, Debug)]
pub(crate) enum Command {
    /// Summarise a structure.
    Info {
        /// The file to read.
        input: PathBuf,
        /// Break the summary down per chain.
        #[arg(long)]
        detail: bool,
    },
    /// Convert a structure from one format to another.
    Convert {
        /// The file to read.
        input: PathBuf,
        /// Where to write.
        output: PathBuf,
        /// Rename a chain on the way out, as `FROM=TO`. May be repeated.
        #[arg(long = "chain-map", value_name = "FROM=TO")]
        chain_map: Vec<String>,
        /// Write counts past the field limits in the extended scheme.
        #[arg(long)]
        hybrid36: bool,
        /// Preserve lossless CIF categories and source ordering.
        #[arg(long)]
        preserve: bool,
        /// Explicit mmCIF/ModelCIF data-block identifier when `_entry.id` is absent.
        #[arg(long, value_name = "ID")]
        cif_block_id: Option<String>,
        /// Allow deterministic connection identifiers for the lowered bond graph.
        #[arg(long)]
        generate_cif_connection_ids: bool,
        /// Explicit mmCIF dictionary connection type for graph bonds.
        #[arg(long, value_name = "TYPE")]
        cif_connection_type: Option<String>,
    },
    /// Check a structure against the invariants it must satisfy.
    Validate {
        /// The file to read.
        input: PathBuf,
        /// Checks to run; multiple names may be comma-separated.
        #[arg(long, value_enum, value_delimiter = ',', default_value = "core")]
        checks: Vec<ValidationChoice>,
        /// CCD required by chemistry-aware geometry checks.
        #[command(flatten)]
        chemistry: CcdArguments,
        /// Maximum bond-length departure in Angstrom.
        #[arg(long)]
        bond_tolerance: Option<f32>,
        /// Maximum aromatic-plane RMS departure in Angstrom.
        #[arg(long)]
        planarity_tolerance: Option<f64>,
        /// Relative convergence tolerance for aromatic plane fitting.
        #[arg(long)]
        plane_relative_tolerance: Option<f64>,
        /// Hard ceiling on cyclic-Jacobi sweeps for aromatic plane fitting.
        #[arg(long)]
        plane_maximum_sweeps: Option<usize>,
        /// Tolerated van der Waals overlap in Angstrom.
        #[arg(long)]
        clash_tolerance: Option<f32>,
        /// Published radius set for clash validation.
        #[arg(long, value_enum)]
        radii: Option<RadiusChoice>,
        /// Expected occupancy sum for each alternate-site atom group.
        #[arg(long)]
        altloc_expected_sum: Option<f64>,
        /// Accepted absolute departure from the expected alternate occupancy sum.
        #[arg(long)]
        altloc_tolerance: Option<f64>,
        /// Absolute z-score threshold for B-factor outlier reporting.
        #[arg(long)]
        b_factor_z_score: Option<f64>,
    },
    /// Measure a structure: extent, radius of gyration, centre.
    Measure {
        /// The file to read.
        input: PathBuf,
    },
    /// Compare two structures atom by atom, after fitting one onto the other.
    Rmsd {
        /// The structure to move.
        mobile: PathBuf,
        /// The structure to move it onto.
        reference: PathBuf,
        /// Report the deviation where the structures sit, without fitting.
        #[arg(long)]
        no_fit: bool,
        /// Declarative atom selection evaluated in both structures.
        #[arg(long)]
        on: Option<String>,
    },
    /// Print the policy an analysis runs under by default.
    Policy,
    /// Inspect built-in analysis profiles.
    Profile {
        #[command(subcommand)]
        command: ProfileCommand,
    },
    /// Intrinsic and backbone geometry workflows.
    Geometry {
        #[command(subcommand)]
        command: GeometryCommand,
    },
    /// Build and characterize a solvent-excluded surface.
    Surface {
        #[command(flatten)]
        args: SurfaceArguments,
    },
    /// Dense multi-model ensemble workflows.
    Ensemble {
        #[command(subcommand)]
        command: EnsembleCommand,
    },
    /// Inspect and convert coordinate trajectories.
    #[command(alias = "traj")]
    Trajectory {
        #[command(subcommand)]
        command: TrajectoryCommand,
    },
    /// Inspect system containers that carry topology and coordinates.
    System {
        #[command(subcommand)]
        command: SystemCommand,
    },
    /// List CIF data blocks and categories with exact counts.
    Categories { input: PathBuf },
    /// Extract, convert, align and compare biological sequences.
    Sequence {
        #[command(subcommand)]
        command: SequenceCommand,
    },
    /// List biological assemblies without materialising coordinates.
    Assemblies { input: PathBuf },
    /// Show rows from one lossless CIF category.
    Show {
        input: PathBuf,
        #[arg(long)]
        category: String,
        #[arg(long)]
        limit: Option<usize>,
    },
    /// List atom contacts at an explicit cutoff in Angstrom.
    Contacts {
        input: PathBuf,
        /// Restrict contacts to two declarative selections.
        #[arg(long, num_args = 2, value_names = ["SELECTION_A", "SELECTION_B"])]
        between: Option<Vec<String>>,
        #[arg(long)]
        cutoff: f32,
    },
    /// List neighbours of a declarative atom selection.
    Neighbors {
        input: PathBuf,
        query: String,
        #[arg(long)]
        cutoff: f32,
    },
    /// Assign secondary structure with the native DSSP kernel.
    Sse {
        input: PathBuf,
        #[command(flatten)]
        chemistry: CcdArguments,
        #[arg(long)]
        electrostatic_prefactor: f64,
        #[arg(long)]
        hydrogen_bond_energy: f64,
        #[arg(long)]
        amide_hydrogen_distance: f32,
        #[arg(long)]
        minimum_sequence_separation: usize,
        #[arg(long)]
        helix_offset: usize,
        #[arg(long, num_args = 2, value_names = ["MIN", "MAX"])]
        turn_offsets: Vec<usize>,
    },
    /// Report backbone torsions with explicit missing values.
    Torsions {
        input: PathBuf,
        #[command(flatten)]
        chemistry: CcdArguments,
    },
    /// List interface residues between two named chains.
    Interfaces {
        input: PathBuf,
        #[arg(long, num_args = 2, value_names = ["CHAIN_A", "CHAIN_B"])]
        between: Vec<String>,
        #[arg(long)]
        cutoff: f32,
    },
    /// Compute solvent-accessible area with explicit physical policy.
    Sasa {
        input: PathBuf,
        #[arg(long)]
        probe: f32,
        #[arg(long)]
        points: u16,
        #[arg(long, value_enum)]
        radii: RadiusChoice,
    },
    /// Compare corresponding coordinates with selected native metrics.
    Compare {
        model: PathBuf,
        reference: PathBuf,
        #[arg(long, value_enum, value_delimiter = ',', num_args = 1..)]
        metrics: Vec<MetricChoice>,
        /// Inclusion radius in Angstrom; required for lDDT.
        #[arg(long)]
        lddt_radius: Option<f32>,
        /// Exclude reference distances at or below this value for lDDT.
        #[arg(long)]
        lddt_minimum_distance: Option<f64>,
        /// Explicit lDDT error tolerances in Angstrom.
        #[arg(long, value_delimiter = ',', num_args = 1..)]
        lddt_tolerances: Vec<f64>,
        /// Empty-domain behavior for lDDT.
        #[arg(long, value_enum)]
        lddt_empty: Option<EmptyLddtChoice>,
        /// Receptor chain used by `DockQ`.
        #[arg(long)]
        receptor: Option<String>,
        /// Ligand chain used by `DockQ`.
        #[arg(long)]
        ligand: Option<String>,
        /// Native-contact distance in Angstrom used by `DockQ`.
        #[arg(long)]
        contact_distance: Option<f32>,
        /// Ligand RMSD scale in Angstrom used by `DockQ`.
        #[arg(long)]
        ligand_scale: Option<f64>,
        /// Interface RMSD scale in Angstrom used by `DockQ`.
        #[arg(long)]
        interface_scale: Option<f64>,
    },
    /// Report semantic structural and metadata differences.
    Diff {
        left: PathBuf,
        right: PathBuf,
        /// Maximum coordinate displacement considered unchanged, in Angstrom.
        #[arg(long)]
        coordinate_tolerance: f32,
    },
    /// Materialise a declarative atom selection or stream its atom rows.
    Select {
        input: PathBuf,
        query: String,
        /// Write the selected structure; its suffix selects the format.
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Print only the selected atom count.
        #[arg(long, conflicts_with = "output")]
        count: bool,
    },
    /// Rigidly align a structure using an explicit selection query.
    Superpose {
        mobile: PathBuf,
        reference: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long)]
        on: String,
        /// Relative threshold used to reject collinear selections.
        #[arg(long)]
        collinear_relative_tolerance: f64,
        /// Relative convergence tolerance for the quaternion eigensolver.
        #[arg(long)]
        eigen_relative_tolerance: f64,
        /// Hard ceiling on complete cyclic-Jacobi sweeps.
        #[arg(long)]
        eigen_maximum_sweeps: usize,
    },
    /// Map polymer chains by sequence with explicit scoring policy.
    MapChains {
        model: PathBuf,
        reference: PathBuf,
        #[command(flatten)]
        chemistry: CcdArguments,
        #[arg(long)]
        min_identity: f64,
        #[arg(long)]
        match_score: i32,
        #[arg(long)]
        mismatch_score: i32,
        #[arg(long)]
        gap_open: i32,
        #[arg(long)]
        gap_extend: i32,
    },
    /// Generate a shell completion script from this command definition.
    Completions { shell: CompletionShell },
    /// Inspect an explicitly versioned local Chemical Component Dictionary.
    Ccd {
        #[command(subcommand)]
        command: CcdCommand,
    },
    /// Retrieve and validate one structure from an explicit endpoint template.
    Fetch {
        id: String,
        /// URL containing the literal placeholder {id}.
        #[arg(long)]
        url_template: String,
        /// Pinned SHA-256 digest of the expected response.
        #[arg(long)]
        sha256: String,
        /// Maximum accepted response size in bytes.
        #[arg(long)]
        max_bytes: u64,
        /// Whole-request timeout in seconds.
        #[arg(long)]
        timeout_seconds: u64,
        /// Maximum HTTP redirects; zero refuses redirects.
        #[arg(long)]
        redirect_limit: usize,
        /// Exact destination file; existing files are not replaced.
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Audit selection stability over explicit policy alternatives.
    Audit {
        #[command(flatten)]
        args: AuditArguments,
    },
    /// Declarative functional-geometry evaluation.
    Fx {
        #[command(subcommand)]
        command: FxCommand,
    },
    /// Deterministic parallel execution over paths and glob patterns.
    Batch {
        #[command(subcommand)]
        command: BatchCommand,
    },
    /// Generate manual pages from the live command definition.
    Man {
        /// Existing directory receiving the generated manual page.
        #[arg(long)]
        outdir: PathBuf,
    },
}

#[derive(Clone, Copy, Subcommand, Debug)]
pub(crate) enum ProfileCommand {
    /// List stable built-in profile identifiers.
    List,
}

#[derive(Subcommand, Debug)]
pub(crate) enum CcdCommand {
    /// Get one component definition.
    Get {
        component: String,
        #[command(flatten)]
        chemistry: CcdArguments,
    },
    /// Retrieve, verify and validate a complete CCD release.
    Update {
        #[arg(long)]
        url: String,
        #[arg(long)]
        sha256: String,
        #[arg(long)]
        max_bytes: u64,
        #[arg(long)]
        timeout_seconds: u64,
        #[arg(long)]
        redirect_limit: usize,
        #[arg(long)]
        version: String,
        #[arg(short, long)]
        output: PathBuf,
        /// Atomically replace an existing destination.
        #[arg(long)]
        replace: bool,
    },
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let context = match cli.global.context() {
        Ok(context) => context,
        Err(error) => {
            eprintln!("policy configuration failed: {error}");
            return exit_code(exit::Exit::Policy);
        }
    };
    if let Err(error) = context.prepare_outputs() {
        eprintln!("output preparation failed: {error}");
        return ExitCode::FAILURE;
    }

    exit_code(dispatch::execute(cli.command, context))
}

fn exit_code(exit: exit::Exit) -> ExitCode {
    // Every documented code fits a byte; anything else would be a defect here.
    match u8::try_from(exit.code()) {
        Ok(code) => ExitCode::from(code),
        Err(_) => ExitCode::FAILURE,
    }
}

#[cfg(test)]
#[path = "main_tests.rs"]
mod tests;

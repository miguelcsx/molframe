//! Arguments for grouped structure, ensemble, trajectory and system workflows.

use crate::RadiusChoice;
use clap::{Args, Subcommand};
use std::path::PathBuf;

/// Optional per-command override for the configured chemistry resource.
#[derive(Args, Debug)]
pub(crate) struct CcdArguments {
    /// Chemical Component Dictionary file; requires --ccd-version.
    #[arg(long, requires = "ccd_version")]
    pub(crate) ccd: Option<PathBuf>,
    /// Exact CCD release identifier; requires --ccd.
    #[arg(long, requires = "ccd")]
    pub(crate) ccd_version: Option<String>,
}

#[derive(Subcommand, Debug)]
pub(crate) enum TrajectoryCommand {
    /// Summarise frames, atoms, time, cells and optional streams.
    Info {
        /// Coordinate trajectory; its suffix selects the reader.
        input: PathBuf,
        /// Optional structure used to verify the atom count.
        #[arg(short, long)]
        topology: Option<PathBuf>,
        /// Angstrom represented by one GSD length unit (required for GSD).
        #[arg(long, value_name = "ANGSTROM_PER_UNIT")]
        gsd_length_scale: Option<f64>,
    },
    /// Convert coordinate trajectories while preserving available frame data.
    Convert {
        /// Coordinate trajectory; its suffix selects the reader.
        input: PathBuf,
        /// Destination trajectory; its suffix selects the writer.
        output: PathBuf,
        /// Keep every Nth frame.
        #[arg(long, default_value_t = 1, value_parser = positive_usize)]
        stride: usize,
        /// Angstrom represented by one GSD length unit (required for GSD).
        #[arg(long, value_name = "ANGSTROM_PER_UNIT")]
        gsd_length_scale: Option<f64>,
        /// TRZ title (required when converting a non-TRZ source to TRZ).
        #[arg(long, value_name = "TEXT")]
        trz_title: Option<String>,
    },
    /// Extract, reorder or repeat selected frames into a trajectory.
    Extract {
        input: PathBuf,
        output: PathBuf,
        #[arg(long, value_delimiter = ',', num_args = 1..)]
        frames: Vec<usize>,
        #[arg(long, value_name = "ANGSTROM_PER_UNIT")]
        gsd_length_scale: Option<f64>,
        #[arg(long, value_name = "TEXT")]
        trz_title: Option<String>,
    },
    /// Measure every frame against one reference frame.
    Rmsd {
        input: PathBuf,
        #[arg(long)]
        reference: usize,
        #[arg(long)]
        no_fit: bool,
        #[arg(long, value_name = "ANGSTROM_PER_UNIT")]
        gsd_length_scale: Option<f64>,
    },
}

fn positive_usize(value: &str) -> Result<usize, String> {
    match value.parse::<usize>() {
        Ok(number) if number > 0 => Ok(number),
        _ => Err("value must be a positive integer".to_owned()),
    }
}

#[derive(Subcommand, Debug)]
pub(crate) enum SystemCommand {
    /// Summarise a DMS system database.
    Info {
        /// DMS system database.
        input: PathBuf,
    },
    /// Copy a DMS system through the validated native model.
    Copy {
        /// DMS system database to read.
        input: PathBuf,
        /// New DMS destination; existing files are never replaced.
        output: PathBuf,
    },
}

#[derive(Subcommand, Debug)]
pub(crate) enum GeometryCommand {
    /// Report backbone torsions with explicit missing values.
    Torsions {
        input: PathBuf,
        #[command(flatten)]
        chemistry: CcdArguments,
    },
    /// Report discrete C-alpha and HELANAL-compatible geometry.
    Helix {
        input: PathBuf,
        #[command(flatten)]
        chemistry: CcdArguments,
        /// Relative convergence tolerance for the symmetric eigensolver.
        #[arg(long)]
        eigen_relative_tolerance: f64,
        /// Hard ceiling on complete cyclic-Jacobi sweeps.
        #[arg(long)]
        eigen_maximum_sweeps: usize,
    },
    /// Build and characterize an indexed solvent-excluded surface.
    Surface {
        #[command(flatten)]
        args: SurfaceArguments,
    },
}

/// Complete explicit definition of a solvent-excluded-surface workflow.
#[derive(Args, Debug)]
pub(crate) struct SurfaceArguments {
    pub(crate) input: PathBuf,
    #[arg(long)]
    pub(crate) probe: f32,
    #[arg(long)]
    pub(crate) resolution: f32,
    /// Hard bound on allocated surface-grid cells.
    #[arg(long)]
    pub(crate) max_cells: usize,
    #[arg(long, value_enum)]
    pub(crate) radii: RadiusChoice,
    #[arg(long)]
    pub(crate) source_vertex: Option<u32>,
    #[arg(long)]
    pub(crate) patch_radius: Option<f64>,
    /// Optional Wavefront OBJ destination for the indexed mesh.
    #[arg(short, long)]
    pub(crate) output: Option<PathBuf>,
}

#[derive(Subcommand, Debug)]
pub(crate) enum EnsembleCommand {
    /// Cartesian PCA over dense models.
    Pca {
        input: PathBuf,
        #[arg(long)]
        components: usize,
        #[arg(long)]
        memory_limit: usize,
        #[arg(long)]
        no_fit: bool,
        /// Iterative fit convergence tolerance; required unless --no-fit.
        #[arg(long)]
        fit_tolerance: Option<f64>,
        /// Iterative fit limit; required unless --no-fit.
        #[arg(long)]
        fit_iterations: Option<usize>,
    },
    /// PCA over backbone torsions using cosine/sine features.
    TorsionPca {
        input: PathBuf,
        #[arg(long)]
        components: usize,
        #[arg(long)]
        memory_limit: usize,
        #[command(flatten)]
        chemistry: CcdArguments,
    },
    /// Diffusion map over fitted-RMSD distances.
    Diffusion {
        input: PathBuf,
        #[arg(long)]
        epsilon: f64,
        #[arg(long)]
        time: u32,
        #[arg(long)]
        dimensions: usize,
        #[arg(long)]
        memory_limit: usize,
    },
}

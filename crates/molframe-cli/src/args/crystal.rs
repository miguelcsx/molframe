//! Arguments for crystallographic data commands.

use clap::{Args, Subcommand, ValueEnum};
use std::path::PathBuf;

/// What `mtz-info` prints.
#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum MtzTable {
    /// Header properties as `property`, `value` rows.
    Summary,
    /// One row per column: label, type and dataset.
    Columns,
}

/// A density histogram: its bin count and the density range it covers.
#[derive(Args, Debug)]
pub(crate) struct HistogramArguments {
    /// Number of equally wide bins; with --histogram-range, adds histogram rows.
    #[arg(long, requires = "histogram_range")]
    pub(crate) histogram_bins: Option<usize>,
    /// Inclusive density range of the histogram.
    #[arg(long, num_args = 2, value_names = ["MIN", "MAX"], allow_hyphen_values = true, requires = "histogram_bins")]
    pub(crate) histogram_range: Option<Vec<f32>>,
}

/// Crystal neighbour search over the structure's recorded symmetry.
#[derive(Args, Debug)]
pub(crate) struct MatesArguments {
    /// A structure whose space group is recorded (`CRYST1` or mmCIF symmetry).
    pub(crate) input: PathBuf,
    /// Largest atom-to-atom distance to a symmetry mate, in Angstrom.
    #[arg(long)]
    pub(crate) cutoff: f64,
    /// Bound on candidate lattice images examined.
    #[arg(long)]
    pub(crate) candidate_limit: usize,
}

#[derive(Subcommand, Debug)]
pub(crate) enum CrystalCommand {
    /// Summarise an MTZ reflection file.
    MtzInfo {
        /// The MTZ file.
        input: PathBuf,
        /// What to print.
        #[arg(long, value_enum)]
        table: MtzTable,
    },
    /// Mean, sigma and extrema of an MRC/CCP4 density map.
    MapStats {
        /// The MRC or CCP4 map.
        input: PathBuf,
        #[command(flatten)]
        histogram: HistogramArguments,
    },
    /// Atoms of neighbouring symmetry copies within a cutoff.
    Mates {
        #[command(flatten)]
        args: MatesArguments,
    },
    /// Pearson correlation of two identically placed density maps.
    MapCorrelation {
        /// The observed map.
        observed: PathBuf,
        /// The calculated map on the same grid.
        calculated: PathBuf,
    },
}

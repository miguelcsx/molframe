//! Arguments for small-molecule chemistry commands.

use super::workflows::CcdArguments;
use clap::{Args, Subcommand, ValueEnum};
use std::path::PathBuf;

/// The named PEOE parameter collection.
#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum PeoeProfileChoice {
    /// Hybridisation-aware Gasteiger-Marsili parameters.
    GasteigerMarsili,
}

impl From<PeoeProfileChoice> for molframe::chemistry::PeoeParameterProfile {
    fn from(value: PeoeProfileChoice) -> Self {
        match value {
            PeoeProfileChoice::GasteigerMarsili => Self::GasteigerMarsili,
        }
    }
}

/// Partial charges by iterative partial equalisation of orbital electronegativity.
#[derive(Args, Debug)]
pub(crate) struct PeoeArguments {
    /// The file to read.
    pub(crate) input: PathBuf,
    #[command(flatten)]
    pub(crate) chemistry: CcdArguments,
    /// Number of charge-transfer iterations.
    #[arg(long)]
    pub(crate) passes: usize,
    /// Damping applied on the first iteration.
    #[arg(long)]
    pub(crate) initial_damping: f64,
    /// Multiplier applied to the damping after every iteration.
    #[arg(long)]
    pub(crate) damping_factor: f64,
    /// Electronegativity differences at or below this transfer no charge.
    #[arg(long)]
    pub(crate) minimum_difference: f64,
    /// Parameter collection.
    #[arg(long, value_enum)]
    pub(crate) profile: PeoeProfileChoice,
}

/// A SMARTS substructure query.
#[derive(Args, Debug)]
pub(crate) struct SmartsArguments {
    /// The structure to search, or with `--component` the component identifier.
    pub(crate) target: String,
    #[command(flatten)]
    pub(crate) chemistry: CcdArguments,
    /// The SMARTS pattern.
    #[arg(long)]
    pub(crate) pattern: String,
    /// Treat the target as a component identifier in the dictionary.
    #[arg(long)]
    pub(crate) component: bool,
}

#[derive(Subcommand, Debug)]
pub(crate) enum ChemCommand {
    /// Gasteiger-Marsili partial charges per atom.
    Peoe {
        #[command(flatten)]
        args: PeoeArguments,
    },
    /// Substructure matches of a SMARTS pattern.
    Smarts {
        #[command(flatten)]
        args: SmartsArguments,
    },
}

//! Controls for scoring structures whose chains and atoms do not line up.

use super::CcdArguments;
use clap::Args;

/// Correspondence controls shared by the scores that can map chains first.
///
/// Every control is required with `--map-chains`, because each decides which
/// atoms are scored.
#[derive(Args, Debug)]
pub(crate) struct MappingArguments {
    /// Match chains by sequence, then residues and atoms, before `DockQ`.
    #[arg(long)]
    pub(crate) map_chains: bool,
    #[command(flatten)]
    pub(crate) chemistry: CcdArguments,
    /// Smallest sequence identity at which two chains correspond, 0 to 1.
    #[arg(long, requires = "map_chains")]
    pub(crate) min_identity: Option<f64>,
    /// Alignment score of identical residues.
    #[arg(long, requires = "map_chains", allow_hyphen_values = true)]
    pub(crate) match_score: Option<i32>,
    /// Alignment score of different residues.
    #[arg(long, requires = "map_chains", allow_hyphen_values = true)]
    pub(crate) mismatch_score: Option<i32>,
    /// Alignment penalty for opening a gap.
    #[arg(long, requires = "map_chains", allow_hyphen_values = true)]
    pub(crate) gap_open: Option<i32>,
    /// Alignment penalty for extending a gap.
    #[arg(long, requires = "map_chains", allow_hyphen_values = true)]
    pub(crate) gap_extend: Option<i32>,
    /// Most chemically equivalent atom mappings tried per residue.
    #[arg(long, requires = "map_chains")]
    pub(crate) automorphism_limit: Option<usize>,
}

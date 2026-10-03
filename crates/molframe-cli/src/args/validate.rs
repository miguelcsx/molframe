//! Arguments for `validate`: which checks run and every threshold they need.
//!
//! A check that embodies a scientific definition takes its thresholds as
//! flags with no default, and refuses to run without them.

use super::workflows::CcdArguments;
use crate::{RadiusChoice, ValidationChoice};
use clap::Args;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub(crate) struct ValidateArguments {
    /// The file to read.
    pub(crate) input: PathBuf,
    /// Checks to run; multiple names may be comma-separated.
    #[arg(long, value_enum, value_delimiter = ',', default_value = "core")]
    pub(crate) checks: Vec<ValidationChoice>,
    /// CCD required by chemistry-aware checks.
    #[command(flatten)]
    pub(crate) chemistry: CcdArguments,
    /// Maximum bond-length departure in Angstrom.
    #[arg(long)]
    pub(crate) bond_tolerance: Option<f32>,
    /// Maximum aromatic-plane RMS departure in Angstrom.
    #[arg(long)]
    pub(crate) planarity_tolerance: Option<f64>,
    /// Relative convergence tolerance for plane fitting.
    #[arg(long)]
    pub(crate) plane_relative_tolerance: Option<f64>,
    /// Hard ceiling on cyclic-Jacobi sweeps for plane fitting.
    #[arg(long)]
    pub(crate) plane_maximum_sweeps: Option<usize>,
    /// Tolerated van der Waals overlap in Angstrom.
    #[arg(long)]
    pub(crate) clash_tolerance: Option<f32>,
    /// Published radius set for clash validation.
    #[arg(long, value_enum)]
    pub(crate) radii: Option<RadiusChoice>,
    /// Expected occupancy sum for each alternate-site atom group.
    #[arg(long)]
    pub(crate) altloc_expected_sum: Option<f64>,
    /// Accepted absolute departure from the expected alternate occupancy sum.
    #[arg(long)]
    pub(crate) altloc_tolerance: Option<f64>,
    /// Absolute z-score threshold for B-factor outlier reporting.
    #[arg(long)]
    pub(crate) b_factor_z_score: Option<f64>,
    #[command(flatten)]
    pub(crate) polymer: PolymerArguments,
    #[command(flatten)]
    pub(crate) references: ReferenceArguments,
}

/// Inputs of the checks that read polymer backbone and nucleotide geometry.
#[derive(Args, Debug)]
pub(crate) struct PolymerArguments {
    /// Caller-authored JSON/TOML polymer roles with `profile_id` and rules.
    #[arg(long)]
    pub(crate) role_profile: Option<PathBuf>,
    /// Largest |omega| counted as a cis peptide bond, in degrees.
    #[arg(long)]
    pub(crate) cis_threshold_degrees: Option<f64>,
    /// Smallest absolute stereocentre volume treated as non-degenerate.
    #[arg(long)]
    pub(crate) chirality_minimum_volume: Option<f64>,
    /// Maximum ligand bond-length departure in Angstrom.
    #[arg(long)]
    pub(crate) ligand_bond_tolerance: Option<f32>,
    /// Largest RMS departure of a nucleotide base from its plane, in Angstrom.
    #[arg(long)]
    pub(crate) base_plane_deviation: Option<f64>,
    /// Accepted glycosidic bond length interval, in Angstrom.
    #[arg(long, num_args = 2, value_names = ["MIN", "MAX"])]
    pub(crate) glycosidic_range: Option<Vec<f64>>,
    /// Accepted phosphodiester link length interval, in Angstrom.
    #[arg(long, num_args = 2, value_names = ["MIN", "MAX"])]
    pub(crate) phosphodiester_range: Option<Vec<f64>>,
}

/// Inputs of the checks that read caller-supplied references and restraints.
#[derive(Args, Debug)]
pub(crate) struct ReferenceArguments {
    /// Reference library (JSON/TOML) holding the empirical distributions.
    #[arg(long)]
    pub(crate) reference_library: Option<PathBuf>,
    /// Ramachandran basin as `REGION=DISTRIBUTION`; REGION is alpha-right,
    /// beta or alpha-left. Repeat for each basin, in tie-breaking order.
    #[arg(long, value_name = "REGION=DISTRIBUTION")]
    pub(crate) basin: Vec<BasinArgument>,
    /// Minimum empirical support of a Ramachandran bin, in [0, 1].
    #[arg(long)]
    pub(crate) ramachandran_minimum_probability: Option<f64>,
    /// Rotamer profile (JSON/TOML) naming each torsion and its distribution.
    #[arg(long)]
    pub(crate) rotamer_profile: Option<PathBuf>,
    /// Minimum empirical support of a rotamer bin, in [0, 1].
    #[arg(long)]
    pub(crate) rotamer_minimum_probability: Option<f64>,
    /// Plane restraints (JSON/TOML) with an id and a selection each.
    #[arg(long)]
    pub(crate) plane_restraints: Option<PathBuf>,
    /// TLS groups (JSON/TOML) with selection, origin and tensors.
    #[arg(long)]
    pub(crate) tls_groups: Option<PathBuf>,
    /// Largest |observed - predicted| B factor accepted under TLS, in square Angstrom.
    #[arg(long)]
    pub(crate) tls_max_deviation: Option<f64>,
    /// Tolerance for the symmetry of the TLS tensors.
    #[arg(long)]
    pub(crate) tls_symmetry_tolerance: Option<f64>,
}

/// One Ramachandran basin and the reference grid that describes it.
#[derive(Clone, Debug)]
pub(crate) struct BasinArgument {
    pub(crate) region: molframe::validation::RamachandranRegion,
    pub(crate) distribution: String,
}

impl std::str::FromStr for BasinArgument {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        use molframe::validation::RamachandranRegion as Region;
        let Some((region, distribution)) = value.split_once('=') else {
            return Err("expected REGION=DISTRIBUTION".to_owned());
        };
        let region = match region {
            "alpha-right" => Region::AlphaHelixRight,
            "beta" => Region::BetaSheet,
            "alpha-left" => Region::AlphaHelixLeft,
            _ => return Err("REGION must be alpha-right, beta or alpha-left".to_owned()),
        };
        if distribution.is_empty() {
            return Err("DISTRIBUTION must not be empty".to_owned());
        }
        Ok(Self {
            region,
            distribution: distribution.to_owned(),
        })
    }
}

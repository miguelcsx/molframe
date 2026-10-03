//! Arguments for the non-covalent interaction analyses.
//!
//! Every geometric threshold is a required flag: the library refuses to carry a
//! literature default, so the command line cannot either.

use super::workflows::CcdArguments;
use clap::{Args, Subcommand};
use std::path::PathBuf;

/// The structure an interaction analysis reads and the dictionary that
/// annotates it.
#[derive(Args, Debug)]
pub(crate) struct InteractionInput {
    /// The file to read.
    pub(crate) input: PathBuf,
    #[command(flatten)]
    pub(crate) chemistry: CcdArguments,
}

/// Geometry of one donor-hydrogen-acceptor hydrogen bond.
#[derive(Args, Debug)]
pub(crate) struct HydrogenBondArguments {
    #[command(flatten)]
    pub(crate) source: InteractionInput,
    /// Largest donor to acceptor heavy-atom distance, in Angstrom.
    #[arg(long)]
    pub(crate) max_donor_acceptor_distance: f32,
    /// Smallest donor-hydrogen-acceptor angle, in degrees.
    #[arg(long)]
    pub(crate) min_angle_degrees: f64,
    /// Apply the structure's unit cell with the minimum-image convention.
    #[arg(long)]
    pub(crate) periodic: bool,
}

/// Distance between oppositely charged atoms.
#[derive(Args, Debug)]
pub(crate) struct SaltBridgeArguments {
    #[command(flatten)]
    pub(crate) source: InteractionInput,
    /// Largest anion to cation distance, in Angstrom.
    #[arg(long)]
    pub(crate) max_distance: f32,
}

/// Numerical controls for the aromatic ring plane fit.
#[derive(Args, Debug)]
pub(crate) struct PlaneFitArguments {
    /// Relative convergence tolerance of the cyclic-Jacobi plane fit.
    #[arg(long)]
    pub(crate) eigen_relative_tolerance: f64,
    /// Ceiling on complete cyclic-Jacobi sweeps of the plane fit.
    #[arg(long)]
    pub(crate) eigen_maximum_sweeps: usize,
}

/// Ring-ring stacking geometry.
#[derive(Args, Debug)]
pub(crate) struct PiStackingArguments {
    #[command(flatten)]
    pub(crate) source: InteractionInput,
    /// Largest distance between ring centres, in Angstrom.
    #[arg(long)]
    pub(crate) max_centre_distance: f32,
    /// Largest acute plane angle classified as parallel, in degrees.
    #[arg(long)]
    pub(crate) max_parallel_angle: f64,
    /// Smallest acute plane angle classified as T-shaped, in degrees.
    #[arg(long)]
    pub(crate) min_perpendicular_angle: f64,
    #[command(flatten)]
    pub(crate) plane_fit: PlaneFitArguments,
}

/// Cation-ring geometry.
#[derive(Args, Debug)]
pub(crate) struct CationPiArguments {
    #[command(flatten)]
    pub(crate) source: InteractionInput,
    /// Largest cation to ring-centre distance, in Angstrom.
    #[arg(long)]
    pub(crate) max_distance: f32,
    /// Largest angle between the ring normal and the cation direction, in degrees.
    #[arg(long)]
    pub(crate) max_face_angle: f64,
    #[command(flatten)]
    pub(crate) plane_fit: PlaneFitArguments,
}

/// Non-covalent interaction analyses over CCD-annotated chemistry.
#[derive(Subcommand, Debug)]
pub(crate) enum InteractionCommand {
    /// Oriented hydrogen bonds from CCD donor and acceptor roles.
    Hbonds {
        #[command(flatten)]
        args: HydrogenBondArguments,
    },
    /// Oppositely charged atom pairs within a distance.
    SaltBridges {
        #[command(flatten)]
        args: SaltBridgeArguments,
    },
    /// Parallel and T-shaped stacking of aromatic rings.
    PiStacking {
        #[command(flatten)]
        args: PiStackingArguments,
    },
    /// Cations over the face of an aromatic ring.
    CationPi {
        #[command(flatten)]
        args: CationPiArguments,
    },
    /// Solvent atoms bridging two partners through hydrogen bonds.
    WaterBridges {
        #[command(flatten)]
        args: HydrogenBondArguments,
    },
}

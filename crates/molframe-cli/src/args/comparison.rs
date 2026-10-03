//! Parameters of the comparison metrics beyond the coordinate-wise scores.
//!
//! A metric that embodies a scientific definition takes every parameter as a
//! flag with no default and refuses to run without it.

use crate::RadiusChoice;
use clap::{Args, ValueEnum};

/// Whether the model is fitted onto the reference before a region is measured.
#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum AlignmentChoice {
    /// Fit the model onto the reference over the measured atoms.
    Fit,
    /// Measure in the frame the inputs are already in.
    None,
}

#[derive(Args, Debug)]
pub(crate) struct ExtraMetricArguments {
    /// Smallest residue-index separation of a counted contact.
    #[arg(long)]
    pub(crate) min_separation: Option<u32>,
    /// Solvent probe radius of the contact-area surface, in Angstrom (CAD).
    #[arg(long)]
    pub(crate) cad_probe: Option<f32>,
    /// Surface point density of the contact-area surface (CAD).
    #[arg(long)]
    pub(crate) cad_density: Option<f32>,
    /// Van der Waals radius set used by the contact-area surface (CAD).
    #[arg(long, value_enum)]
    pub(crate) radii: Option<RadiusChoice>,
    /// Selection of the guide atoms CE aligns, evaluated in each structure.
    #[arg(long)]
    pub(crate) guide: Option<String>,
    /// Consecutive guide atoms per aligned CE fragment pair.
    #[arg(long)]
    pub(crate) ce_window: Option<usize>,
    /// Largest insertion CE examines between fragments.
    #[arg(long)]
    pub(crate) ce_max_gap: Option<usize>,
    /// Number of candidate CE paths retained.
    #[arg(long)]
    pub(crate) ce_max_paths: Option<usize>,
    /// Fragment-pair similarity bound CE accepts; it must exceed the path bound.
    #[arg(long, allow_hyphen_values = true)]
    pub(crate) ce_fragment_threshold: Option<f64>,
    /// Path similarity bound CE accepts.
    #[arg(long, allow_hyphen_values = true)]
    pub(crate) ce_path_threshold: Option<f64>,
    /// Byte ceiling of the CE search workspace.
    #[arg(long)]
    pub(crate) ce_memory_limit: Option<usize>,
    /// Report the original CE significance (z-score) calibration.
    #[arg(long)]
    pub(crate) ce_significance: bool,
    /// Selection of the interface or pocket atoms, evaluated in the reference.
    #[arg(long)]
    pub(crate) region: Option<String>,
    /// Whether the model is fitted onto the reference over the region.
    #[arg(long, value_enum)]
    pub(crate) alignment: Option<AlignmentChoice>,
    /// Selection of the ligand atoms, evaluated in the reference.
    #[arg(long)]
    pub(crate) ligand_selection: Option<String>,
    /// CCD component identifier of the ligand.
    #[arg(long)]
    pub(crate) ligand_component: Option<String>,
    /// Most chemically equivalent ligand atom mappings tried.
    #[arg(long)]
    pub(crate) ligand_automorphism_limit: Option<usize>,
    /// Optional selection the model is fitted onto before the ligand is measured.
    #[arg(long)]
    pub(crate) align_on: Option<String>,
}

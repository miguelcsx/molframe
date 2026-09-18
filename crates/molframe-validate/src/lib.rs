//! Geometric and completeness validation.
//!
//! Each check reads a structure and returns a sorted list of the specific atoms
//! or residues that fail it, not just a count — a validation that only says how
//! many problems exist cannot be acted on. The findings are deterministic and
//! ordered on stable keys so two runs agree exactly.

#![forbid(unsafe_code)]

#[path = "planarity.rs"]
mod aromatic_planarity;
mod backbone;
#[path = "ramachandran.rs"]
mod backbone_conformation;
#[path = "completeness.rs"]
mod chain_integrity;
#[path = "ccd_completeness.rs"]
mod component_completeness;
#[path = "occupancy.rs"]
mod conformer_occupancy;
#[path = "valence.rs"]
mod coordination_valence;
mod covalent;
#[path = "bond_length.rs"]
mod covalent_length;
#[path = "reference.rs"]
mod distributions;
#[path = "plane_restraint.rs"]
mod explicit_planes;
#[path = "reference_geometry.rs"]
mod geometric_reference;
#[path = "ligand.rs"]
mod ligand_geometry;
#[path = "real_space.rs"]
mod map_correlation;
#[path = "nucleic.rs"]
mod nucleic_geometry;
mod numeric;
#[path = "peptide.rs"]
mod peptide_geometry;
#[path = "governed/mod.rs"]
mod policy_execution;
#[path = "quality.rs"]
mod scalar_quality;
#[path = "rotamer.rs"]
mod sidechain_conformation;
#[path = "chirality.rs"]
mod stereocentre;
#[path = "clash.rs"]
mod steric_overlap;
#[path = "bfactor/mod.rs"]
mod thermal_motion;

pub use aromatic_planarity::{
    PlanarityError, PlanarityFlag, PlanarityOptions, nonplanar_aromatic_rings,
};
pub use backbone_conformation::{
    RamachandranBasin, RamachandranError, RamachandranOptions, RamachandranRecord,
    RamachandranRegion, classify, ramachandran, ramachandran_outliers,
};
pub use chain_integrity::{ChainCompleteness, CompletenessError, completeness};
pub use component_completeness::{
    CcdCompletenessReport, ResidueAtomCompleteness, ccd_missing_atoms,
};
pub use conformer_occupancy::{
    AltlocOccupancyError, AltlocOccupancyIssue, AltlocOccupancyOptions, AltlocOccupancyRecord,
    AltlocOccupancyReport, altloc_occupancy_sums,
};
pub use coordination_valence::{ValenceError, overvalent_atoms};
pub use covalent_length::{BondDeviation, bond_length_deviations};
pub use distributions::{
    ReferenceAssessment, ReferenceDistribution, ReferenceError, ReferenceLibrary,
    assess_bond_deviation, assess_ramachandran,
};
pub use explicit_planes::{
    PlaneRestraint, PlaneRestraintFlag, PlaneRestraintReport, plane_restraint_outliers,
};
pub use geometric_reference::{
    ReferenceAngleFlag, ReferenceBondFlag, ReferenceGeometryOptions, ReferenceGeometryReport,
    reference_geometry,
};
pub use ligand_geometry::{LigandGeometryReport, ligand_geometry, ligand_geometry_outliers};
pub use map_correlation::{
    RealSpaceCorrelation, RealSpaceCorrelationError, masked_real_space_correlation,
    real_space_map_correlation, sampled_real_space_correlation,
};
pub use nucleic_geometry::{
    NucleicGeometryError, NucleicGeometryIssue, NucleicGeometryPolicy, NucleicGeometryRecord,
    nucleic_acid_geometry,
};
pub use peptide_geometry::{CisPeptide, cis_peptides};
pub use policy_execution::{
    AltlocOccupancyKernelError, BFactorKernelError, CcdCompletenessKernelError, GovernedMapError,
    LigandGeometryKernelError, PlaneRestraintKernelError, altloc_occupancy_sums_kernel,
    b_factor_distribution_kernel, bond_length_deviations_kernel, ccd_missing_atoms_kernel,
    chirality_kernel, cis_peptides_kernel, clashes_kernel, completeness_kernel,
    governed_masked_real_space_correlation, governed_real_space_map_correlation,
    governed_sampled_real_space_correlation, ligand_geometry_kernel, nucleic_geometry_kernel,
    planarity_kernel, plane_restraint_outliers_kernel, quality_flags_kernel, ramachandran_kernel,
    ramachandran_outliers_kernel, reference_geometry_kernel, rotamer_kernel,
    tls_b_factor_consistency_kernel, valence_kernel,
};
pub use scalar_quality::{QualityFlag, QualityIssue, quality_flags};
pub use sidechain_conformation::{
    RotamerDefinition, RotamerError, RotamerFlag, RotamerOptions, RotamerProfile, RotamerReport,
    rotamer_outliers,
};
pub use stereocentre::{
    ChiralityFlag, ChiralityIssue, ChiralityOptions, ChiralityReport, chirality_outliers,
};
pub use steric_overlap::{Clash, clashes};
pub use thermal_motion::{
    BFactorDistribution, BFactorError, BFactorOutlier, TlsBFactorFlag, TlsBFactorReport, TlsGroup,
    TlsModel, b_factor_distribution, tls_b_factor_consistency,
};

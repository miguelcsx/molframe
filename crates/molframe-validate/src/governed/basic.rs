use crate::{
    AltlocOccupancyError, AltlocOccupancyOptions, AltlocOccupancyReport, BondDeviation,
    CcdCompletenessReport, ChainCompleteness, CisPeptide, Clash, CompletenessError,
    LigandGeometryReport, PlanarityError, PlanarityFlag, PlanarityOptions, QualityFlag,
    ValenceError, altloc_occupancy_sums, bond_length_deviations, ccd_missing_atoms, cis_peptides,
    clashes, completeness, ligand_geometry, nonplanar_aromatic_rings, overvalent_atoms,
    quality_flags,
};
use molframe_analysis::{AnalysisDescriptor, FrameKernelResult, StructureKernel, structure_kernel};
use molframe_chem::RadiusSet;
use molframe_core::contract::{AnalysisPolicy, Coverage, ParameterValue, Status};
use molframe_core::{Diagnostic, ExecutionContext, Structure};
use molframe_spatial::{SpatialBackend, SpatialError};
use std::convert::Infallible;

fn descriptor(name: &'static str) -> AnalysisDescriptor {
    AnalysisDescriptor::new(name, "1")
}

fn float(value: impl Into<f64>) -> ParameterValue {
    ParameterValue::Float(value.into().to_bits())
}

fn integer(value: usize) -> ParameterValue {
    match i64::try_from(value) {
        Ok(value) => ParameterValue::Integer(value),
        Err(_) => ParameterValue::Text(value.to_string().into()),
    }
}

fn complete<T>(structure: &Structure, value: T) -> FrameKernelResult<T> {
    FrameKernelResult::complete(value, structure.atom_count())
}

/// Governed ligand-geometry setup failure.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum LigandGeometryKernelError {
    /// Selected bond counts exceed public coverage counters.
    #[error("ligand bond count exceeds u32 coverage limits")]
    CoverageOverflow,
}

/// Governed alternate-occupancy setup failure.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum AltlocOccupancyKernelError {
    /// The raw validation kernel refused its input.
    #[error(transparent)]
    Validation(#[from] AltlocOccupancyError),
    /// Coverage counts exceed the public contract.
    #[error("alternate occupancy count exceeds u32 coverage limits")]
    CoverageOverflow,
}

/// Governed CCD completeness setup failure.
#[derive(Debug, thiserror::Error)]
pub enum CcdCompletenessKernelError {
    /// CCD or identifier failure.
    #[error("CCD completeness diagnostic: {0:?}")]
    Diagnostic(Diagnostic),
    /// Coverage counts exceed the public contract.
    #[error("CCD completeness count exceeds u32 coverage limits")]
    CoverageOverflow,
}

impl From<Diagnostic> for CcdCompletenessKernelError {
    fn from(value: Diagnostic) -> Self {
        Self::Diagnostic(value)
    }
}

/// Governed steric-clash validation.
#[must_use]
pub fn clashes_kernel(
    tolerance: f32,
    radii: RadiusSet,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = Vec<Clash>, Error = SpatialError> {
    structure_kernel(
        descriptor("steric-clashes")
            .with_parameter("tolerance", float(tolerance))
            .with_parameter(
                "radius_set",
                ParameterValue::Text(format!("{radii:?}").into()),
            )
            .with_parameter(
                "spatial_backend",
                ParameterValue::Text(format!("{spatial:?}").into()),
            ),
        move |structure: &Structure, _policy: &AnalysisPolicy, context: &ExecutionContext| {
            clashes(structure, tolerance, radii, spatial, context)
                .map(|value| complete(structure, value))
        },
    )
}

/// Governed bond-length validation.
#[must_use]
pub fn bond_length_deviations_kernel(
    tolerance: f32,
) -> impl StructureKernel<Output = Vec<BondDeviation>, Error = Infallible> {
    structure_kernel(
        descriptor("bond-length-deviations").with_parameter("tolerance", float(tolerance)),
        move |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            Ok(complete(
                structure,
                bond_length_deviations(structure, tolerance),
            ))
        },
    )
}

/// Governed hetero-component bond-geometry validation.
#[must_use]
pub fn ligand_geometry_kernel(
    tolerance: f32,
) -> impl StructureKernel<Output = LigandGeometryReport, Error = LigandGeometryKernelError> {
    structure_kernel(
        descriptor("ligand-geometry").with_parameter("tolerance", float(tolerance)),
        move |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            let value = ligand_geometry(structure, tolerance);
            let intended = u32::try_from(value.intended)
                .map_err(|_| LigandGeometryKernelError::CoverageOverflow)?;
            let used = u32::try_from(value.assessed)
                .map_err(|_| LigandGeometryKernelError::CoverageOverflow)?;
            let missing = intended
                .checked_sub(used)
                .ok_or(LigandGeometryKernelError::CoverageOverflow)?;
            Ok(FrameKernelResult::governed(
                value,
                if missing == 0 {
                    Status::Complete
                } else {
                    Status::Partial
                },
                Coverage {
                    intended,
                    used,
                    missing,
                    ambiguous: 0,
                },
            ))
        },
    )
}

/// Governed cis-peptide validation.
#[must_use]
pub fn cis_peptides_kernel(
    threshold_degrees: f64,
) -> impl StructureKernel<Output = Vec<CisPeptide>, Error = Diagnostic> {
    structure_kernel(
        descriptor("cis-peptides").with_parameter("threshold_degrees", float(threshold_degrees)),
        move |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            cis_peptides(structure, threshold_degrees).map(|value| complete(structure, value))
        },
    )
}

/// Governed aromatic-planarity validation.
#[must_use]
pub fn planarity_kernel(
    options: PlanarityOptions,
) -> impl StructureKernel<Output = Vec<PlanarityFlag>, Error = PlanarityError> {
    structure_kernel(
        descriptor("aromatic-planarity")
            .with_parameter("maximum_deviation", float(options.maximum_deviation))
            .with_parameter(
                "plane_fit_relative_tolerance",
                float(options.plane_fit.relative_tolerance),
            )
            .with_parameter(
                "plane_fit_maximum_sweeps",
                integer(options.plane_fit.maximum_sweeps),
            ),
        move |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            nonplanar_aromatic_rings(structure, options).map(|value| complete(structure, value))
        },
    )
}

/// Governed occupancy and displacement-factor validation.
#[must_use]
pub fn quality_flags_kernel() -> impl StructureKernel<Output = Vec<QualityFlag>, Error = Infallible>
{
    structure_kernel(
        descriptor("coordinate-quality-flags"),
        move |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            Ok(complete(structure, quality_flags(structure)))
        },
    )
}

/// Governed alternate-location occupancy-sum validation.
#[must_use]
pub fn altloc_occupancy_sums_kernel(
    options: AltlocOccupancyOptions,
) -> impl StructureKernel<Output = AltlocOccupancyReport, Error = AltlocOccupancyKernelError> {
    structure_kernel(
        descriptor("altloc-occupancy-sums")
            .with_parameter("expected_sum", float(options.expected_sum))
            .with_parameter("tolerance", float(options.tolerance)),
        move |structure: &Structure, policy: &AnalysisPolicy, _context: &ExecutionContext| {
            let value = altloc_occupancy_sums(structure, policy.identifiers, options)?;
            let intended = u32::try_from(value.intended)
                .map_err(|_| AltlocOccupancyKernelError::CoverageOverflow)?;
            let used = u32::try_from(value.assessed)
                .map_err(|_| AltlocOccupancyKernelError::CoverageOverflow)?;
            let missing = intended
                .checked_sub(used)
                .ok_or(AltlocOccupancyKernelError::CoverageOverflow)?;
            Ok(FrameKernelResult::governed(
                value,
                if missing == 0 {
                    Status::Complete
                } else {
                    Status::Partial
                },
                Coverage {
                    intended,
                    used,
                    missing,
                    ambiguous: 0,
                },
            ))
        },
    )
}

/// Governed per-residue CCD missing-atom validation.
#[must_use]
pub fn ccd_missing_atoms_kernel(
    provider: &dyn molframe_chem::ComponentProvider,
) -> impl StructureKernel<Output = CcdCompletenessReport, Error = CcdCompletenessKernelError> + '_ {
    structure_kernel(
        descriptor("ccd-missing-atoms"),
        move |structure: &Structure, policy: &AnalysisPolicy, _context: &ExecutionContext| {
            let value = ccd_missing_atoms(structure, provider, policy)?;
            let intended = u32::try_from(value.intended)
                .map_err(|_| CcdCompletenessKernelError::CoverageOverflow)?;
            let used = u32::try_from(value.assessed)
                .map_err(|_| CcdCompletenessKernelError::CoverageOverflow)?;
            let ambiguous = u32::try_from(value.ambiguous)
                .map_err(|_| CcdCompletenessKernelError::CoverageOverflow)?;
            let missing = intended
                .checked_sub(used)
                .and_then(|count| count.checked_sub(ambiguous))
                .ok_or(CcdCompletenessKernelError::CoverageOverflow)?;
            Ok(FrameKernelResult::governed(
                value,
                if missing == 0 && ambiguous == 0 {
                    Status::Complete
                } else {
                    Status::Partial
                },
                Coverage {
                    intended,
                    used,
                    missing,
                    ambiguous,
                },
            ))
        },
    )
}

/// Governed chain-completeness validation.
#[must_use]
pub fn completeness_kernel()
-> impl StructureKernel<Output = Vec<ChainCompleteness>, Error = CompletenessError> {
    structure_kernel(
        descriptor("chain-completeness"),
        move |structure: &Structure, policy: &AnalysisPolicy, _context: &ExecutionContext| {
            completeness(structure, policy.identifiers).map(|value| complete(structure, value))
        },
    )
}

/// Governed valence validation.
#[must_use]
pub fn valence_kernel() -> impl StructureKernel<Output = Vec<ValenceError>, Error = Infallible> {
    structure_kernel(
        descriptor("overvalent-atoms"),
        move |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            Ok(complete(structure, overvalent_atoms(structure)))
        },
    )
}

use crate::{
    ChiralityOptions, ChiralityReport, NucleicGeometryError, NucleicGeometryPolicy,
    NucleicGeometryRecord, RamachandranError, RamachandranOptions, RamachandranRecord,
    ReferenceGeometryOptions, ReferenceGeometryReport, ReferenceLibrary, RotamerError,
    RotamerOptions, RotamerProfile, RotamerReport, chirality_outliers, nucleic_acid_geometry,
    ramachandran, reference_geometry, rotamer_outliers,
};
use molframe_analysis::{AnalysisDescriptor, FrameKernelResult, StructureKernel, structure_kernel};
use molframe_chem::{ComponentProvider, PolymerRoleProfile};
use molframe_core::contract::{AltlocPolicy, AnalysisPolicy, Coverage, ParameterValue, Status};
use molframe_core::{Code, Diagnostic, ExecutionContext, Structure};

fn float(value: f64) -> ParameterValue {
    ParameterValue::Float(value.to_bits())
}

fn integer(value: usize) -> ParameterValue {
    match i64::try_from(value) {
        Ok(value) => ParameterValue::Integer(value),
        Err(_) => ParameterValue::Text(value.to_string().into()),
    }
}

/// Governed CCD-driven nucleic-geometry validation.
#[must_use]
pub fn nucleic_geometry_kernel<'a>(
    provider: &'a dyn ComponentProvider,
    roles: &'a PolymerRoleProfile,
    options: NucleicGeometryPolicy,
) -> impl StructureKernel<Output = Vec<NucleicGeometryRecord>, Error = NucleicGeometryError> + 'a {
    structure_kernel(
        AnalysisDescriptor::new("nucleic-acid-geometry", "1")
            .with_parameter(
                "maximum_base_plane_deviation",
                float(options.maximum_base_plane_deviation),
            )
            .with_parameter(
                "glycosidic_bond_minimum",
                float(options.glycosidic_bond_range[0]),
            )
            .with_parameter(
                "glycosidic_bond_maximum",
                float(options.glycosidic_bond_range[1]),
            )
            .with_parameter(
                "phosphodiester_bond_minimum",
                float(options.phosphodiester_bond_range[0]),
            )
            .with_parameter(
                "phosphodiester_bond_maximum",
                float(options.phosphodiester_bond_range[1]),
            )
            .with_parameter(
                "plane_fit_relative_tolerance",
                float(options.plane_fit.relative_tolerance),
            )
            .with_parameter(
                "plane_fit_maximum_sweeps",
                integer(options.plane_fit.maximum_sweeps),
            ),
        move |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            nucleic_acid_geometry(structure, provider, roles, options)
                .map(|value| FrameKernelResult::complete(value, structure.atom_count()))
        },
    )
}

fn resolved_policy(policy: &AnalysisPolicy) -> AnalysisPolicy {
    policy.clone().with_altloc(AltlocPolicy::KeepAll)
}

fn chirality_result(
    value: ChiralityReport,
) -> Result<FrameKernelResult<ChiralityReport>, Diagnostic> {
    let intended = u32::try_from(value.intended)
        .map_err(|_| Diagnostic::new(Code::E1901).with_context("limit", "coverage counters"))?;
    let used = u32::try_from(value.assessed)
        .map_err(|_| Diagnostic::new(Code::E1901).with_context("limit", "coverage counters"))?;
    let missing = intended
        .checked_sub(used)
        .ok_or_else(|| Diagnostic::new(Code::E3011).with_context("field", "chirality coverage"))?;
    let mut result = FrameKernelResult::governed(
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
    );
    result.warnings.clone_from(&result.value.findings);
    Ok(result)
}

fn ramachandran_result(
    structure: &Structure,
    value: Vec<RamachandranRecord>,
) -> Result<FrameKernelResult<Vec<RamachandranRecord>>, RamachandranError> {
    let intended = structure
        .data()
        .chains()
        .map(|chain| {
            let residues = chain.residues().count();
            match residues.checked_sub(2) {
                Some(interior) => interior,
                None => 0,
            }
        })
        .sum::<usize>();
    let intended = u32::try_from(intended).map_err(|_| RamachandranError::CoverageOverflow)?;
    let used = u32::try_from(value.len()).map_err(|_| RamachandranError::CoverageOverflow)?;
    let missing = intended
        .checked_sub(used)
        .ok_or(RamachandranError::CoverageOverflow)?;
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
}

fn rotamer_result(value: RotamerReport) -> Result<FrameKernelResult<RotamerReport>, RotamerError> {
    let intended = u32::try_from(value.intended).map_err(|_| RotamerError::CoverageOverflow)?;
    let used = u32::try_from(value.assessed).map_err(|_| RotamerError::CoverageOverflow)?;
    let missing = intended
        .checked_sub(used)
        .ok_or(RotamerError::CoverageOverflow)?;
    let mut result = FrameKernelResult::governed(
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
    );
    result.warnings.clone_from(&result.value.findings);
    Ok(result)
}

fn reference_geometry_result(
    value: ReferenceGeometryReport,
) -> Result<FrameKernelResult<ReferenceGeometryReport>, Diagnostic> {
    let intended = u32::try_from(value.intended)
        .map_err(|_| Diagnostic::new(Code::E1901).with_context("limit", "coverage counters"))?;
    let used = u32::try_from(value.assessed)
        .map_err(|_| Diagnostic::new(Code::E1901).with_context("limit", "coverage counters"))?;
    let missing = intended.checked_sub(used).ok_or_else(|| {
        Diagnostic::new(Code::E3011).with_context("field", "reference geometry coverage")
    })?;
    let mut result = FrameKernelResult::governed(
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
    );
    result.warnings.clone_from(&result.value.findings);
    Ok(result)
}

/// Governed CCD-declared stereocentre validation.
#[must_use]
pub fn chirality_kernel(
    provider: &dyn ComponentProvider,
    options: ChiralityOptions,
) -> impl StructureKernel<Output = ChiralityReport, Error = molframe_core::Diagnostic> + '_ {
    structure_kernel(
        AnalysisDescriptor::new("chirality", "1")
            .with_parameter("minimum_abs_volume", float(options.minimum_abs_volume)),
        move |structure: &Structure, policy: &AnalysisPolicy, _context: &ExecutionContext| {
            chirality_outliers(structure, provider, &resolved_policy(policy), options)
                .and_then(chirality_result)
        },
    )
}

/// Governed CCD ideal-coordinate bond and angle validation.
#[must_use]
pub fn reference_geometry_kernel(
    provider: &dyn ComponentProvider,
    options: ReferenceGeometryOptions,
) -> impl StructureKernel<Output = ReferenceGeometryReport, Error = Diagnostic> + '_ {
    structure_kernel(
        AnalysisDescriptor::new("reference-geometry", "1")
            .with_parameter(
                "maximum_bond_deviation",
                float(options.maximum_bond_deviation),
            )
            .with_parameter(
                "maximum_angle_deviation_degrees",
                float(options.maximum_angle_deviation_degrees),
            ),
        move |structure: &Structure, policy: &AnalysisPolicy, _context: &ExecutionContext| {
            reference_geometry(structure, provider, policy.identifiers, options)
                .and_then(reference_geometry_result)
        },
    )
}

/// Governed full Ramachandran assessment against explicit reference grids.
#[must_use]
pub fn ramachandran_kernel(
    options: RamachandranOptions<'_>,
) -> impl StructureKernel<Output = Vec<RamachandranRecord>, Error = RamachandranError> + '_ {
    structure_kernel(
        AnalysisDescriptor::new("ramachandran", "1")
            .with_parameter("minimum_probability", float(options.minimum_probability()))
            .with_parameter(
                "reference_id",
                ParameterValue::Text(options.references().id().into()),
            )
            .with_parameter(
                "reference_version",
                ParameterValue::Text(options.references().version().into()),
            ),
        move |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            ramachandran(structure, &options)
                .and_then(|value| ramachandran_result(structure, value))
        },
    )
}

/// Governed Ramachandran outlier assessment.
#[must_use]
pub fn ramachandran_outliers_kernel(
    options: RamachandranOptions<'_>,
) -> impl StructureKernel<Output = Vec<RamachandranRecord>, Error = RamachandranError> + '_ {
    structure_kernel(
        AnalysisDescriptor::new("ramachandran-outliers", "1")
            .with_parameter("minimum_probability", float(options.minimum_probability()))
            .with_parameter(
                "reference_id",
                ParameterValue::Text(options.references().id().into()),
            )
            .with_parameter(
                "reference_version",
                ParameterValue::Text(options.references().version().into()),
            ),
        move |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            let assessed = ramachandran(structure, &options)?;
            let outliers = assessed
                .iter()
                .filter(|record| record.region == crate::RamachandranRegion::Outlier)
                .cloned()
                .collect();
            let mut result = ramachandran_result(structure, assessed)?;
            result.value = outliers;
            Ok(result)
        },
    )
}

/// Governed CCD-connected rotamer validation against an explicit reference library.
#[must_use]
pub fn rotamer_kernel<'a>(
    provider: &'a dyn ComponentProvider,
    references: &'a ReferenceLibrary,
    profile: &'a RotamerProfile,
    options: RotamerOptions,
) -> impl StructureKernel<Output = RotamerReport, Error = RotamerError> + 'a {
    structure_kernel(
        AnalysisDescriptor::new("rotamer-outliers", "1")
            .with_parameter("minimum_probability", float(options.minimum_probability))
            .with_parameter("reference_id", ParameterValue::Text(references.id().into()))
            .with_parameter(
                "reference_version",
                ParameterValue::Text(references.version().into()),
            )
            .with_parameter("profile_id", ParameterValue::Text(profile.id().into()))
            .with_parameter(
                "profile_version",
                ParameterValue::Text(profile.version().into()),
            ),
        move |structure: &Structure, policy: &AnalysisPolicy, _context: &ExecutionContext| {
            rotamer_outliers(
                structure,
                provider,
                &resolved_policy(policy),
                references,
                profile,
                options,
            )
            .and_then(rotamer_result)
        },
    )
}

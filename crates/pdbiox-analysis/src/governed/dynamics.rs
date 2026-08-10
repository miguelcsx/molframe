//! Governed entry points for trajectory-derived dynamics observables.

use super::common::{descriptor, float, integer};
use super::standalone::StandaloneAnalysisError;
use crate::{
    DielectricOptions, DielectricResult, DynamicsError, WaterDynamicsOptions, WaterLag,
    dielectric_from_dipoles, water_dynamics,
};
use pdbiox_core::AtomSelection;
use pdbiox_core::contract::{
    Analysis, AnalysisPolicy, Coverage, ParameterValue, Provenance, SourceRef, Status,
};

fn complete<T>(
    value: T,
    observations: usize,
    policy: &AnalysisPolicy,
    descriptor: &super::AnalysisDescriptor,
) -> Result<Analysis<T>, StandaloneAnalysisError<DynamicsError>> {
    let observations =
        u32::try_from(observations).map_err(|_| StandaloneAnalysisError::CoverageOverflow)?;
    Ok(Analysis {
        value,
        status: Status::Complete,
        coverage: Coverage::complete(observations),
        warnings: Vec::new(),
        assumptions: Vec::new(),
        provenance: descriptor.apply(Provenance::new(policy).with_source(SourceRef::Memory)),
    })
}

/// Computes water survival and residence with explicit stable water IDs.
///
/// # Errors
///
/// Returns a dynamics-domain error or coverage overflow.
pub fn governed_water_dynamics(
    occupancy: &[AtomSelection],
    options: WaterDynamicsOptions,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<WaterLag>>, StandaloneAnalysisError<DynamicsError>> {
    let value = water_dynamics(occupancy, options).map_err(StandaloneAnalysisError::Kernel)?;
    complete(
        value,
        occupancy.len(),
        policy,
        &descriptor("water-dynamics")
            .with_parameter("maximum_lag", integer(options.maximum_lag))
            .with_parameter("identity", ParameterValue::Text("caller-stable-ids".into())),
    )
}

/// Computes dielectric response with an explicit unit/boundary prefactor.
///
/// # Errors
///
/// Returns a dynamics-domain error or coverage overflow.
pub fn governed_dielectric_from_dipoles(
    dipoles: &[[f64; 3]],
    options: DielectricOptions,
    policy: &AnalysisPolicy,
) -> Result<Analysis<DielectricResult>, StandaloneAnalysisError<DynamicsError>> {
    let value =
        dielectric_from_dipoles(dipoles, options).map_err(StandaloneAnalysisError::Kernel)?;
    complete(
        value,
        dipoles.len(),
        policy,
        &descriptor("dielectric-fluctuation")
            .with_parameter("volume", float(options.volume))
            .with_parameter("temperature", float(options.temperature))
            .with_parameter(
                "fluctuation_prefactor",
                float(options.fluctuation_prefactor),
            ),
    )
}

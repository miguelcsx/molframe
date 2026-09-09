//! Policy-bound physical analysis operations.

use super::requests::{ExecutionPlanError, PlanInput};
use pdbiox_core::contract::{Analysis, AnalysisPolicy, PeriodicPolicy};
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;
use pdbiox_spatial::SpatialBackend;

/// A reusable structure-bound physical analysis request.
#[derive(Clone, Debug)]
pub enum PhysicalRequest {
    /// Normalized site-site radial distribution.
    RadialDistribution {
        /// Sites contributing the first side of the distribution.
        left: AtomSelection,
        /// Sites contributing the second side of the distribution.
        right: AtomSelection,
        /// Shell and normalization controls.
        options: pdbiox_analysis::RadialDistributionOptions,
        /// Compatibility switch for minimum-image periodicity.
        periodic: bool,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Neighbour counts for every selected left site.
    CoordinationNumbers {
        /// Sites whose neighbours are counted.
        left: AtomSelection,
        /// Sites accepted as neighbours.
        right: AtomSelection,
        /// Inclusive lower shell bound.
        minimum_distance: f32,
        /// Exclusive upper shell bound.
        maximum_distance: f32,
        /// Spatial implementation.
        backend: SpatialBackend,
        /// Compatibility switch for minimum-image periodicity.
        periodic: bool,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Connected components of explicit membrane representative sites.
    Leaflets {
        /// Representative sites supplied by the caller.
        sites: AtomSelection,
        /// Component construction controls.
        options: pdbiox_analysis::LeafletOptions,
        /// Compatibility switch for minimum-image periodicity.
        periodic: bool,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Weighted linear density along one Cartesian axis.
    LinearDensity {
        /// Slot containing atom-aligned weights.
        weights: usize,
        /// Axis, bounds and binning controls.
        options: pdbiox_analysis::LinearDensityOptions,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Weighted Cartesian density grid.
    DensityMap {
        /// Slot containing atom-aligned weights.
        weights: usize,
        /// Grid origin, spacing and shape.
        spec: pdbiox_analysis::DensityGridSpec,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Pore profile over explicit atom radii.
    PoreProfile {
        /// Slot containing atom-aligned radii.
        radii: usize,
        /// Axis and sampling controls.
        options: pdbiox_analysis::PoreProfileOptions,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Exposed-surface contacts over explicit atom radii.
    SurfaceContacts {
        /// Slot containing atom-aligned radii.
        radii: usize,
        /// Surface overlap tolerance.
        tolerance: f32,
        /// Probe radius.
        probe: f32,
        /// Surface sampling density.
        density: f32,
        /// Minimum exposed contact area.
        minimum_area: f32,
        /// Spatial implementation.
        backend: SpatialBackend,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
}

/// A typed result from a [`PhysicalRequest`].
#[derive(Clone, Debug)]
pub enum PhysicalValue {
    /// Radial bins and their analysis contract.
    RadialDistribution(Analysis<Vec<pdbiox_analysis::RadialBin>>),
    /// Per-left-site coordination counts and their analysis contract.
    CoordinationNumbers(Analysis<Vec<u32>>),
    /// Deterministic connected components and their analysis contract.
    Leaflets(Analysis<Vec<pdbiox_analysis::Leaflet>>),
    /// Linear density bins and their analysis contract.
    LinearDensity(Analysis<Vec<pdbiox_analysis::LinearDensityBin>>),
    /// Density grid and its analysis contract.
    DensityMap(Analysis<pdbiox_analysis::DensityGrid>),
    /// Pore samples and their analysis contract.
    PoreProfile(Analysis<Vec<pdbiox_analysis::PoreSample>>),
    /// Exposed-surface contacts and their analysis contract.
    SurfaceContacts(Analysis<Vec<pdbiox_analysis::Contact>>),
}

/// Execute one physical request through the governed analysis layer.
pub fn execute(
    operation: &str,
    request: &PhysicalRequest,
    structure: &Structure,
    input: PlanInput<'_>,
    context: &pdbiox_core::ExecutionContext,
) -> Result<PhysicalValue, ExecutionPlanError> {
    match request {
        PhysicalRequest::RadialDistribution {
            left,
            right,
            options,
            periodic,
            policy,
        } => {
            let policy = effective_policy(policy, *periodic);
            let kernel = pdbiox_analysis::radial_distribution_kernel(left, right, *options);
            pdbiox_analysis::analyse_structure(structure, &policy, &kernel, context)
                .map(PhysicalValue::RadialDistribution)
                .map_err(governed_error)
        }
        PhysicalRequest::CoordinationNumbers {
            left,
            right,
            minimum_distance,
            maximum_distance,
            backend,
            periodic,
            policy,
        } => {
            let policy = effective_policy(policy, *periodic);
            let kernel = pdbiox_analysis::coordination_numbers_kernel(
                left,
                right,
                *minimum_distance,
                *maximum_distance,
                *backend,
            );
            pdbiox_analysis::analyse_structure(structure, &policy, &kernel, context)
                .map(PhysicalValue::CoordinationNumbers)
                .map_err(governed_error)
        }
        PhysicalRequest::Leaflets {
            sites,
            options,
            periodic,
            policy,
        } => {
            let policy = effective_policy(policy, *periodic);
            let kernel = pdbiox_analysis::leaflets_kernel(sites, *options);
            pdbiox_analysis::analyse_structure(structure, &policy, &kernel, context)
                .map(PhysicalValue::Leaflets)
                .map_err(governed_error)
        }
        PhysicalRequest::LinearDensity {
            weights,
            options,
            policy,
        } => execute_linear_density(
            operation, structure, input, *weights, *options, policy, context,
        ),
        PhysicalRequest::DensityMap {
            weights,
            spec,
            policy,
        } => execute_density_map(
            operation, structure, input, *weights, *spec, policy, context,
        ),
        PhysicalRequest::PoreProfile {
            radii,
            options,
            policy,
        } => execute_pore_profile(
            operation, structure, input, *radii, *options, policy, context,
        ),
        PhysicalRequest::SurfaceContacts {
            radii,
            tolerance,
            probe,
            density,
            minimum_area,
            backend,
            policy,
        } => execute_surface_contacts(
            operation,
            structure,
            input,
            *radii,
            SurfaceContactParameters {
                tolerance: *tolerance,
                probe: *probe,
                density: *density,
                minimum_area: *minimum_area,
                backend: *backend,
            },
            policy,
            context,
        ),
    }
}

fn execute_linear_density(
    operation: &str,
    structure: &Structure,
    input: PlanInput<'_>,
    slot: usize,
    options: pdbiox_analysis::LinearDensityOptions,
    policy: &AnalysisPolicy,
    context: &pdbiox_core::ExecutionContext,
) -> Result<PhysicalValue, ExecutionPlanError> {
    let weights = scalar_slot(operation, input, slot)?;
    let kernel = pdbiox_analysis::linear_density_kernel(weights.values, options);
    pdbiox_analysis::analyse_structure(structure, policy, &kernel, context)
        .map(PhysicalValue::LinearDensity)
        .map_err(governed_error)
}

fn execute_density_map(
    operation: &str,
    structure: &Structure,
    input: PlanInput<'_>,
    slot: usize,
    spec: pdbiox_analysis::DensityGridSpec,
    policy: &AnalysisPolicy,
    context: &pdbiox_core::ExecutionContext,
) -> Result<PhysicalValue, ExecutionPlanError> {
    let weights = scalar_slot(operation, input, slot)?;
    let kernel = pdbiox_analysis::density_map_kernel(weights.values, spec);
    pdbiox_analysis::analyse_structure(structure, policy, &kernel, context)
        .map(PhysicalValue::DensityMap)
        .map_err(governed_error)
}

fn execute_pore_profile(
    operation: &str,
    structure: &Structure,
    input: PlanInput<'_>,
    slot: usize,
    options: pdbiox_analysis::PoreProfileOptions,
    policy: &AnalysisPolicy,
    context: &pdbiox_core::ExecutionContext,
) -> Result<PhysicalValue, ExecutionPlanError> {
    let radii = float_slot(operation, input, slot)?;
    let kernel = pdbiox_analysis::pore_profile_kernel(radii.values, options);
    pdbiox_analysis::analyse_structure(structure, policy, &kernel, context)
        .map(PhysicalValue::PoreProfile)
        .map_err(governed_error)
}

#[derive(Clone, Copy)]
struct SurfaceContactParameters {
    tolerance: f32,
    probe: f32,
    density: f32,
    minimum_area: f32,
    backend: SpatialBackend,
}

fn execute_surface_contacts(
    operation: &str,
    structure: &Structure,
    input: PlanInput<'_>,
    slot: usize,
    parameters: SurfaceContactParameters,
    policy: &AnalysisPolicy,
    context: &pdbiox_core::ExecutionContext,
) -> Result<PhysicalValue, ExecutionPlanError> {
    let radii = float_slot(operation, input, slot)?;
    let kernel = pdbiox_analysis::surface_contacts_kernel(
        radii.values,
        pdbiox_analysis::SurfaceContactOptions {
            tolerance: parameters.tolerance,
            probe: parameters.probe,
            surface_density: parameters.density,
            minimum_area: parameters.minimum_area,
            backend: parameters.backend,
        },
    );
    pdbiox_analysis::analyse_structure(structure, policy, &kernel, context)
        .map(PhysicalValue::SurfaceContacts)
        .map_err(governed_error)
}

fn effective_policy(policy: &AnalysisPolicy, periodic: bool) -> AnalysisPolicy {
    let mut effective = policy.clone();
    if periodic {
        effective.periodic = PeriodicPolicy::MinimumImage;
    }
    effective
}

fn governed_error<E: std::fmt::Display>(error: E) -> ExecutionPlanError {
    ExecutionPlanError::Governed(error.to_string().into())
}

fn scalar_slot<'a>(
    operation: &str,
    input: PlanInput<'a>,
    slot: usize,
) -> Result<super::requests::ScalarInput<'a>, ExecutionPlanError> {
    input
        .scalars
        .get(slot)
        .copied()
        .ok_or_else(|| ExecutionPlanError::ScalarSlot {
            operation: operation.into(),
            slot,
        })
}

fn float_slot<'a>(
    operation: &str,
    input: PlanInput<'a>,
    slot: usize,
) -> Result<super::requests::FloatInput<'a>, ExecutionPlanError> {
    input
        .floats
        .get(slot)
        .copied()
        .ok_or_else(|| ExecutionPlanError::FloatSlot {
            operation: operation.into(),
            slot,
        })
}

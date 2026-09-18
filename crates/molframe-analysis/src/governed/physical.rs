//! Policy-aware adapters for atom-aligned physical analyses.

use super::common::{backend, complete, descriptor, float, integer};
use super::{StructureKernel, mapped_structure_kernel};
use crate::{
    CentreGroup, DensityError, DensityGrid, DensityGridSpec, Leaflet, LeafletOptions,
    LinearDensityBin, LinearDensityOptions, PoreError, PoreProfileOptions, PoreSample, RadialBin,
    RadialDistributionOptions, RadialError, centre_of_mass_radial_distribution,
    coordination_numbers, density_map, identify_leaflets, linear_density, pore_profile,
    radial_distribution,
};
use molframe_core::contract::{AnalysisPolicy, ParameterValue, PeriodicPolicy};
use molframe_core::{AtomSelection, ExecutionContext, Structure};
use molframe_spatial::{PeriodicBox, SpatialBackend, SpatialError};

/// Policy projection or scientific-kernel failure.
#[derive(Debug, thiserror::Error)]
pub enum PhysicalKernelError<E> {
    /// The scientific kernel refused its explicit inputs.
    #[error(transparent)]
    Kernel(E),
    /// Periodic analysis was requested without a unit cell.
    #[error("periodic analysis requires a unit cell")]
    MissingCell,
    /// The kernel cannot represent all periodic images distinctly.
    #[error("this kernel supports minimum-image periodicity, not all-image periodicity")]
    AllImagesUnsupported,
    /// An atom-aligned side input does not cover the source topology.
    #[error("atom-aligned side input is shorter than the source topology")]
    SideInputLength,
}

fn periodic<E>(
    structure: &Structure,
    policy: &AnalysisPolicy,
) -> Result<Option<PeriodicBox>, PhysicalKernelError<E>> {
    match policy.periodic {
        PeriodicPolicy::None => Ok(None),
        PeriodicPolicy::MinimumImage => structure
            .data()
            .cell
            .ok_or(PhysicalKernelError::MissingCell)
            .and_then(|cell| {
                PeriodicBox::from_cell(cell)
                    .map(Some)
                    .map_err(|_| PhysicalKernelError::MissingCell)
            }),
        _ => Err(PhysicalKernelError::AllImagesUnsupported),
    }
}

fn selection(source: &AtomSelection, source_atoms: &[usize]) -> AtomSelection {
    AtomSelection::from_sorted(
        source_atoms
            .iter()
            .enumerate()
            .filter_map(|(target, source_atom)| {
                let source_atom = u32::try_from(*source_atom).ok()?;
                source
                    .contains(source_atom)
                    .then(|| u32::try_from(target).ok())?
            })
            .collect(),
    )
}

fn aligned<T: Copy>(values: &[T], source_atoms: &[usize]) -> Option<Vec<T>> {
    source_atoms
        .iter()
        .map(|source| values.get(*source).copied())
        .collect()
}

/// Governed weighted linear-density kernel.
#[must_use]
pub fn linear_density_kernel(
    weights: &[f64],
    options: LinearDensityOptions,
) -> impl StructureKernel<Output = Vec<LinearDensityBin>, Error = PhysicalKernelError<DensityError>> + '_
{
    mapped_structure_kernel(
        descriptor("linear-density")
            .with_parameter(
                "axis",
                ParameterValue::Text(format!("{:?}", options.axis).into()),
            )
            .with_parameter("minimum", float(options.minimum))
            .with_parameter("maximum", float(options.maximum))
            .with_parameter("bins", integer(options.bins)),
        move |structure: &Structure,
              _policy: &AnalysisPolicy,
              source_atoms: &[usize],
              _context: &ExecutionContext| {
            let weights =
                aligned(weights, source_atoms).ok_or(PhysicalKernelError::SideInputLength)?;
            linear_density(structure.positions(), &weights, options)
                .map(|value| complete(structure, value))
                .map_err(PhysicalKernelError::Kernel)
        },
    )
}

/// Governed weighted Cartesian density-map kernel.
#[must_use]
pub fn density_map_kernel(
    weights: &[f64],
    spec: DensityGridSpec,
) -> impl StructureKernel<Output = DensityGrid, Error = PhysicalKernelError<DensityError>> + '_ {
    mapped_structure_kernel(
        descriptor("density-map")
            .with_parameter(
                "origin",
                ParameterValue::Text(format!("{:?}", spec.origin).into()),
            )
            .with_parameter(
                "spacing",
                ParameterValue::Text(format!("{:?}", spec.spacing).into()),
            )
            .with_parameter(
                "shape",
                ParameterValue::Text(format!("{:?}", spec.shape).into()),
            ),
        move |structure: &Structure,
              _policy: &AnalysisPolicy,
              source_atoms: &[usize],
              _context: &ExecutionContext| {
            let weights =
                aligned(weights, source_atoms).ok_or(PhysicalKernelError::SideInputLength)?;
            density_map(structure.positions(), &weights, spec)
                .map(|value| complete(structure, value))
                .map_err(PhysicalKernelError::Kernel)
        },
    )
}

/// Governed membrane-leaflet kernel with policy-driven periodicity.
#[must_use]
pub fn leaflets_kernel(
    sites: &AtomSelection,
    options: LeafletOptions,
) -> impl StructureKernel<Output = Vec<Leaflet>, Error = PhysicalKernelError<SpatialError>> + '_ {
    mapped_structure_kernel(
        descriptor("membrane-leaflets")
            .with_parameter("connection_distance", float(options.connection_distance))
            .with_parameter("spatial_backend", backend(options.backend)),
        move |structure: &Structure,
              policy: &AnalysisPolicy,
              source_atoms: &[usize],
              context: &ExecutionContext| {
            let sites = selection(sites, source_atoms);
            let periodic = periodic(structure, policy)?;
            identify_leaflets(
                structure.positions(),
                &sites,
                options,
                periodic.as_ref(),
                context,
            )
            .map(|value| complete(structure, value))
            .map_err(PhysicalKernelError::Kernel)
        },
    )
}

/// Governed pore-profile kernel over explicit atom radii.
#[must_use]
pub fn pore_profile_kernel(
    radii: &[f32],
    options: PoreProfileOptions,
) -> impl StructureKernel<Output = Vec<PoreSample>, Error = PhysicalKernelError<PoreError>> + '_ {
    mapped_structure_kernel(
        descriptor("pore-profile")
            .with_parameter("samples", integer(options.samples))
            .with_parameter("search_radius", float(options.search_radius))
            .with_parameter("grid_spacing", float(options.grid_spacing))
            .with_parameter("probe_radius", float(options.probe_radius))
            .with_parameter("memory_limit_bytes", integer(options.memory_limit_bytes)),
        move |structure: &Structure,
              _policy: &AnalysisPolicy,
              source_atoms: &[usize],
              _context: &ExecutionContext| {
            let radii = aligned(radii, source_atoms).ok_or(PhysicalKernelError::SideInputLength)?;
            pore_profile(structure.positions(), &radii, options)
                .map(|value| complete(structure, value))
                .map_err(PhysicalKernelError::Kernel)
        },
    )
}

/// Governed radial-distribution kernel.
#[must_use]
pub fn radial_distribution_kernel<'a>(
    left: &'a AtomSelection,
    right: &'a AtomSelection,
    options: RadialDistributionOptions,
) -> impl StructureKernel<Output = Vec<RadialBin>, Error = PhysicalKernelError<RadialError>> + 'a {
    mapped_structure_kernel(
        descriptor("radial-distribution")
            .with_parameter("minimum_distance", float(options.minimum_distance))
            .with_parameter("maximum_distance", float(options.maximum_distance))
            .with_parameter("bins", integer(options.bins))
            .with_parameter("volume", float(options.volume))
            .with_parameter("spatial_backend", backend(options.backend)),
        move |structure: &Structure,
              policy: &AnalysisPolicy,
              source_atoms: &[usize],
              context: &ExecutionContext| {
            let left = selection(left, source_atoms);
            let right = selection(right, source_atoms);
            let periodic = periodic(structure, policy)?;
            radial_distribution(
                structure.positions(),
                &left,
                &right,
                options,
                periodic.as_ref(),
                context,
            )
            .map(|value| complete(structure, value))
            .map_err(PhysicalKernelError::Kernel)
        },
    )
}

/// Governed centre-of-mass radial-distribution kernel.
#[must_use]
pub fn centre_of_mass_radial_distribution_kernel<'a>(
    masses: &'a [f64],
    left: &'a [CentreGroup],
    right: &'a [CentreGroup],
    options: RadialDistributionOptions,
) -> impl StructureKernel<Output = Vec<RadialBin>, Error = PhysicalKernelError<RadialError>> + 'a {
    mapped_structure_kernel(
        descriptor("centre-of-mass-radial-distribution")
            .with_parameter("minimum_distance", float(options.minimum_distance))
            .with_parameter("maximum_distance", float(options.maximum_distance))
            .with_parameter("bins", integer(options.bins))
            .with_parameter("volume", float(options.volume))
            .with_parameter("spatial_backend", backend(options.backend))
            .with_parameter(
                "group_imaging",
                ParameterValue::Text("coordinates-already-whole".into()),
            ),
        move |structure: &Structure,
              policy: &AnalysisPolicy,
              source_atoms: &[usize],
              context: &ExecutionContext| {
            let masses =
                aligned(masses, source_atoms).ok_or(PhysicalKernelError::SideInputLength)?;
            let left = left
                .iter()
                .map(|group| CentreGroup {
                    atoms: selection(&group.atoms, source_atoms),
                })
                .collect::<Vec<_>>();
            let right = right
                .iter()
                .map(|group| CentreGroup {
                    atoms: selection(&group.atoms, source_atoms),
                })
                .collect::<Vec<_>>();
            let periodic = periodic(structure, policy)?;
            centre_of_mass_radial_distribution(
                structure.positions(),
                &masses,
                &left,
                &right,
                options,
                periodic.as_ref(),
                context,
            )
            .map(|value| complete(structure, value))
            .map_err(PhysicalKernelError::Kernel)
        },
    )
}

/// Governed coordination-number kernel.
#[must_use]
pub fn coordination_numbers_kernel<'a>(
    left: &'a AtomSelection,
    right: &'a AtomSelection,
    minimum_distance: f32,
    maximum_distance: f32,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = Vec<u32>, Error = PhysicalKernelError<RadialError>> + 'a {
    mapped_structure_kernel(
        descriptor("coordination-numbers")
            .with_parameter("minimum_distance", float(minimum_distance))
            .with_parameter("maximum_distance", float(maximum_distance))
            .with_parameter("spatial_backend", backend(spatial)),
        move |structure: &Structure,
              policy: &AnalysisPolicy,
              source_atoms: &[usize],
              context: &ExecutionContext| {
            let left = selection(left, source_atoms);
            let right = selection(right, source_atoms);
            let periodic = periodic(structure, policy)?;
            coordination_numbers(
                structure.positions(),
                &left,
                &right,
                crate::CoordinationOptions {
                    minimum_distance,
                    maximum_distance,
                    backend: spatial,
                },
                periodic.as_ref(),
                context,
            )
            .map(|value| complete(structure, value))
            .map_err(PhysicalKernelError::Kernel)
        },
    )
}

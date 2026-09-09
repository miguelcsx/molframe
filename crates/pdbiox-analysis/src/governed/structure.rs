//! Ready-to-run governed kernels for structure-level descriptors.

use super::common::{backend, complete, descriptor, float, integer};
use super::{StructureKernel, structure_kernel};
use crate::{
    DsspError, DsspOptions, GaussianNetworkModel, GnmError, GnmOptions, HalfSphereExposure,
    HseError, NucleicTorsionError, NucleicTorsions, SseRecord, chain_interface,
    gaussian_network_model, half_sphere_exposure, nucleic_torsions, secondary_structure,
};
use pdbiox_core::ExecutionContext;
use pdbiox_core::contract::{AnalysisPolicy, ParameterValue};
use pdbiox_core::index::ResidueIndex;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;
use pdbiox_spatial::{PeriodicBox, SpatialBackend, SpatialError};

/// Governed chain-interface kernel.
#[must_use]
pub fn chain_interface_kernel<'a>(
    first_chain: &'a str,
    second_chain: &'a str,
    cutoff: f32,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = Vec<ResidueIndex>, Error = SpatialError> + 'a {
    structure_kernel(
        descriptor("chain-interface")
            .with_parameter("first_chain", ParameterValue::Text(first_chain.into()))
            .with_parameter("second_chain", ParameterValue::Text(second_chain.into()))
            .with_parameter("cutoff", float(cutoff))
            .with_parameter("spatial_backend", backend(spatial)),
        move |structure: &Structure, _policy: &AnalysisPolicy, context: &ExecutionContext| {
            chain_interface(
                structure,
                first_chain,
                second_chain,
                cutoff,
                spatial,
                context,
            )
            .map(|value| complete(structure, value))
        },
    )
}

/// Governed DSSP-compatible secondary-structure kernel.
#[must_use]
pub fn secondary_structure_kernel(
    options: &DsspOptions,
) -> impl StructureKernel<Output = Vec<SseRecord>, Error = DsspError> + '_ {
    structure_kernel(
        descriptor("secondary-structure-dssp")
            .with_parameter(
                "electrostatic_prefactor",
                float(options.electrostatic_prefactor),
            )
            .with_parameter("hydrogen_bond_energy", float(options.hydrogen_bond_energy))
            .with_parameter(
                "amide_hydrogen_distance",
                float(options.amide_hydrogen_distance),
            )
            .with_parameter(
                "minimum_sequence_separation",
                integer(options.minimum_sequence_separation),
            )
            .with_parameter("helix_offset", integer(options.helix_offset))
            .with_parameter("turn_offset_start", integer(*options.turn_offsets.start()))
            .with_parameter("turn_offset_end", integer(*options.turn_offsets.end())),
        move |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            secondary_structure(structure, options).map(|value| complete(structure, value))
        },
    )
}

/// Governed half-sphere-exposure kernel.
#[must_use]
pub fn half_sphere_exposure_kernel(
    radius: f32,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = Vec<HalfSphereExposure>, Error = HseError> {
    structure_kernel(
        descriptor("half-sphere-exposure")
            .with_parameter("radius", float(radius))
            .with_parameter("spatial_backend", backend(spatial)),
        move |structure: &Structure, _policy: &AnalysisPolicy, context: &ExecutionContext| {
            half_sphere_exposure(structure, radius, spatial, context)
                .map(|value| complete(structure, value))
        },
    )
}

/// Governed nucleic-acid torsion kernel.
#[must_use]
pub fn nucleic_torsions_kernel()
-> impl StructureKernel<Output = Vec<NucleicTorsions>, Error = NucleicTorsionError> {
    structure_kernel(
        descriptor("nucleic-acid-torsions"),
        move |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            nucleic_torsions(structure).map(|value| complete(structure, value))
        },
    )
}

/// Governed Gaussian-network kernel over an explicit site selection.
#[must_use]
pub fn gnm_kernel<'a>(
    sites: &'a AtomSelection,
    options: GnmOptions,
    periodic: Option<&'a PeriodicBox>,
) -> impl StructureKernel<Output = GaussianNetworkModel, Error = GnmError> + 'a {
    structure_kernel(
        descriptor("gaussian-network-model")
            .with_parameter(
                "site_count",
                ParameterValue::Text(sites.len().to_string().into()),
            )
            .with_parameter("contact_distance", float(options.contact_distance))
            .with_parameter("mode_count", integer(options.mode_count))
            .with_parameter("zero_mode_tolerance", float(options.zero_mode_tolerance))
            .with_parameter("memory_limit_bytes", integer(options.memory_limit_bytes))
            .with_parameter("spatial_backend", backend(options.backend))
            .with_parameter("periodic", ParameterValue::Boolean(periodic.is_some())),
        move |structure: &Structure, _policy: &AnalysisPolicy, context: &ExecutionContext| {
            gaussian_network_model(structure.positions(), sites, options, periodic, context)
                .map(|value| complete(structure, value))
        },
    )
}

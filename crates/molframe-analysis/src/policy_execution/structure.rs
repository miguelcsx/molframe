//! Ready-to-run governed kernels for structure-level descriptors.

use super::descriptor_parameters::{backend, complete, descriptor, float, integer};
use super::kernel::FrameKernelResult;
use super::{StructureKernel, structure_kernel};
use crate::{
    DsspError, DsspOptions, GaussianNetworkModel, GnmError, GnmOptions, HalfSphereExposure,
    HseError, NucleicTorsionError, NucleicTorsions, SseTable, chain_interface,
    gaussian_network_model, half_sphere_exposure, nucleic_torsions, protein_domain,
    secondary_structure,
};
use molframe_core::ExecutionContext;
use molframe_core::contract::{AnalysisPolicy, Coverage, ParameterValue, Quality};
use molframe_core::index::ResidueIndex;
use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;
use molframe_spatial::{PeriodicBox, SpatialBackend};

/// Governed chain-interface kernel.
#[must_use]
pub fn chain_interface_kernel<'a>(
    first_chain: &'a str,
    second_chain: &'a str,
    cutoff: f32,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = Vec<ResidueIndex>, Error = crate::InterfaceError> + 'a {
    structure_kernel(
        descriptor("chain-interface")
            .estimating("the residues of one chain that touch another chain within the cutoff")
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
) -> impl StructureKernel<Output = SseTable, Error = DsspError> + '_ {
    structure_kernel(
        descriptor("secondary-structure-dssp")
            .estimating("the DSSP secondary-structure class of each residue from its backbone hydrogen bonds")
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
///
/// A residue needs both an alpha and a beta carbon role — the beta carbon defines
/// the outward direction — so a residue without one takes no part. Those residues
/// are reported as missing coverage rather than silently dropped, which is what
/// lets `missing_atoms` decide whether a partial answer is acceptable: an input
/// whose profile supplies only backbone roles would otherwise return an empty
/// table with no indication that nothing was measured.
#[must_use]
pub fn half_sphere_exposure_kernel(
    radius: f32,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = Vec<HalfSphereExposure>, Error = HseError> {
    structure_kernel(
        descriptor("half-sphere-exposure")
            .estimating("the neighbours of each residue in the half-sphere towards and away from its side chain")
            .requiring(super::Requirement::PolymerSideChainDirection)
            .with_parameter("radius", float(radius))
            .with_parameter("spatial_backend", backend(spatial)),
        move |structure: &Structure, _policy: &AnalysisPolicy, context: &ExecutionContext| {
            let measured = half_sphere_exposure(structure, radius, spatial, context)?;
            // The analysis describes protein burial, so its domain is the
            // residues the profile identified as protein. Counting a bound water
            // or a ligand as an intended protein residue would report an analysis
            // that covered every measurable residue as mostly incomplete, and the
            // caller has no way to tell that from a genuinely truncated run.
            let intended = protein_domain(structure);
            // A structure cannot hold more residues than `u32`, so both counts fit.
            let measured_count = match u32::try_from(measured.len()) {
                Ok(count) => count,
                Err(_) => u32::MAX,
            };
            let used = measured_count.min(intended);
            Ok(FrameKernelResult::governed(
                measured,
                if used < intended {
                    Quality::Partial
                } else {
                    Quality::Complete
                },
                Coverage {
                    intended,
                    used,
                    missing: intended - used,
                    ambiguous: 0,
                },
            ))
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

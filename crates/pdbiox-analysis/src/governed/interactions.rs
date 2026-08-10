//! Ready-to-run governed kernels for intermolecular interactions.

use super::common::{backend, complete, descriptor, float, integer};
use super::{StructureKernel, structure_kernel};
use crate::{
    BasePair, BasePairError, BasePairOptions, CationPi, CationPiError, CationPiOptions, Contact,
    ContactMap, HydrogenBond, HydrogenBondError, HydrogenBondOptions, PiStacking, PiStackingError,
    PiStackingOptions, SaltBridge, WaterBridge, WaterBridgeOptions, atom_contacts, base_pairs,
    cation_pi, hydrogen_bonds, pi_stacking, residue_contact_map, salt_bridges, surface_contacts,
    water_bridges,
};
use pdbiox_chem::ComponentProvider;
use pdbiox_core::contract::{AnalysisPolicy, ParameterValue};
use pdbiox_core::structure::Structure;
use pdbiox_spatial::{SpatialBackend, SpatialError};
use pdbiox_surface::SasaError;

fn with_hydrogen_bond_parameters(
    descriptor: super::AnalysisDescriptor,
    options: HydrogenBondOptions,
) -> super::AnalysisDescriptor {
    descriptor
        .with_parameter(
            "maximum_donor_acceptor_distance",
            float(options.maximum_donor_acceptor_distance),
        )
        .with_parameter(
            "minimum_angle_degrees",
            float(options.minimum_angle_degrees),
        )
        .with_parameter("spatial_backend", backend(options.backend))
        .with_parameter("periodic", ParameterValue::Boolean(options.periodic))
}

fn with_plane_fit(
    descriptor: super::AnalysisDescriptor,
    options: pdbiox_geom::EigenOptions,
) -> super::AnalysisDescriptor {
    descriptor
        .with_parameter(
            "plane_fit_relative_tolerance",
            float(options.relative_tolerance),
        )
        .with_parameter("plane_fit_maximum_sweeps", integer(options.maximum_sweeps))
}

/// Governed atom-contact kernel with complete parameter provenance.
#[must_use]
pub fn contacts_kernel(
    cutoff: f32,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = Vec<Contact>, Error = SpatialError> {
    structure_kernel(
        descriptor("atom-contacts")
            .with_parameter("cutoff", float(cutoff))
            .with_parameter("spatial_backend", backend(spatial)),
        move |structure: &Structure, _policy: &AnalysisPolicy| {
            atom_contacts(structure, cutoff, spatial).map(|value| complete(structure, value))
        },
    )
}

/// Governed residue-contact-map kernel.
#[must_use]
pub fn contact_map_kernel(
    cutoff: f32,
    minimum_separation: u32,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = ContactMap, Error = SpatialError> {
    structure_kernel(
        descriptor("residue-contact-map")
            .with_parameter("cutoff", float(cutoff))
            .with_parameter(
                "minimum_separation",
                ParameterValue::Integer(i64::from(minimum_separation)),
            )
            .with_parameter("spatial_backend", backend(spatial)),
        move |structure: &Structure, _policy: &AnalysisPolicy| {
            residue_contact_map(structure, cutoff, minimum_separation, spatial)
                .map(|value| complete(structure, value))
        },
    )
}

/// Governed CCD-driven hydrogen-bond kernel.
#[must_use]
pub fn hydrogen_bonds_kernel(
    options: HydrogenBondOptions,
) -> impl StructureKernel<Output = Vec<HydrogenBond>, Error = HydrogenBondError> {
    structure_kernel(
        with_hydrogen_bond_parameters(descriptor("hydrogen-bonds"), options),
        move |structure: &Structure, _policy: &AnalysisPolicy| {
            hydrogen_bonds(structure, options).map(|value| complete(structure, value))
        },
    )
}

/// Governed salt-bridge kernel using CCD formal charges.
#[must_use]
pub fn salt_bridges_kernel(
    maximum_distance: f32,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = Vec<SaltBridge>, Error = SpatialError> {
    structure_kernel(
        descriptor("salt-bridges")
            .with_parameter("maximum_distance", float(maximum_distance))
            .with_parameter("spatial_backend", backend(spatial)),
        move |structure: &Structure, _policy: &AnalysisPolicy| {
            salt_bridges(structure, maximum_distance, spatial)
                .map(|value| complete(structure, value))
        },
    )
}

/// Governed aromatic stacking kernel.
#[must_use]
pub fn pi_stacking_kernel(
    options: PiStackingOptions,
) -> impl StructureKernel<Output = Vec<PiStacking>, Error = PiStackingError> {
    structure_kernel(
        with_plane_fit(
            descriptor("pi-stacking")
                .with_parameter(
                    "maximum_centre_distance",
                    float(options.maximum_centre_distance),
                )
                .with_parameter(
                    "maximum_parallel_angle",
                    float(options.maximum_parallel_angle),
                )
                .with_parameter(
                    "minimum_perpendicular_angle",
                    float(options.minimum_perpendicular_angle),
                ),
            options.plane_fit,
        ),
        move |structure: &Structure, _policy: &AnalysisPolicy| {
            pi_stacking(structure, options).map(|value| complete(structure, value))
        },
    )
}

/// Governed CCD cation-pi kernel.
#[must_use]
pub fn cation_pi_kernel(
    options: CationPiOptions,
) -> impl StructureKernel<Output = Vec<CationPi>, Error = CationPiError> {
    structure_kernel(
        with_plane_fit(
            descriptor("cation-pi")
                .with_parameter("maximum_distance", float(options.maximum_distance))
                .with_parameter("maximum_face_angle", float(options.maximum_face_angle)),
            options.plane_fit,
        ),
        move |structure: &Structure, _policy: &AnalysisPolicy| {
            cation_pi(structure, options).map(|value| complete(structure, value))
        },
    )
}

/// Governed water-bridge kernel.
#[must_use]
pub fn water_bridges_kernel(
    options: WaterBridgeOptions,
) -> impl StructureKernel<Output = Vec<WaterBridge>, Error = HydrogenBondError> {
    structure_kernel(
        with_hydrogen_bond_parameters(descriptor("water-bridges"), options.hydrogen_bonds),
        move |structure: &Structure, _policy: &AnalysisPolicy| {
            water_bridges(structure, options).map(|value| complete(structure, value))
        },
    )
}

/// Governed CCD base-pair kernel.
#[must_use]
pub fn base_pairs_kernel(
    provider: &dyn ComponentProvider,
    options: BasePairOptions,
) -> impl StructureKernel<Output = Vec<BasePair>, Error = BasePairError> + '_ {
    structure_kernel(
        with_hydrogen_bond_parameters(
            descriptor("canonical-base-pairs").with_parameter(
                "minimum_hydrogen_bonds",
                integer(options.minimum_hydrogen_bonds),
            ),
            options.hydrogen_bonds,
        ),
        move |structure: &Structure, _policy: &AnalysisPolicy| {
            base_pairs(structure, provider, options).map(|value| complete(structure, value))
        },
    )
}

/// Governed exposed-surface contact kernel.
#[must_use]
pub fn surface_contacts_kernel(
    radii: &[f32],
    tolerance: f32,
    probe: f32,
    density: f32,
    minimum_area: f32,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = Vec<Contact>, Error = SasaError> + '_ {
    structure_kernel(
        descriptor("surface-contacts")
            .with_parameter("tolerance", float(tolerance))
            .with_parameter("probe", float(probe))
            .with_parameter("surface_density", float(density))
            .with_parameter("minimum_area", float(minimum_area))
            .with_parameter("spatial_backend", backend(spatial)),
        move |structure: &Structure, _policy: &AnalysisPolicy| {
            surface_contacts(
                structure,
                radii,
                tolerance,
                probe,
                density,
                minimum_area,
                spatial,
            )
            .map(|value| complete(structure, value))
        },
    )
}

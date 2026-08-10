//! Governed adapters from structure kernels to deterministic frame execution.

mod common;
mod descriptor;
mod dynamics;
mod error;
mod execute;
mod helical;
mod interactions;
mod kernel;
mod native;
mod physical;
mod standalone;
mod structure;

pub use descriptor::AnalysisDescriptor;
pub use dynamics::{governed_dielectric_from_dipoles, governed_water_dynamics};
pub use error::GovernedAnalysisError;
pub use execute::{GovernedStructureAnalysis, analyse_structure, analyse_trajectory};
pub use helical::{governed_helical_parameters, governed_helical_steps};
pub use interactions::{
    base_pairs_kernel, cation_pi_kernel, contact_map_kernel, contacts_kernel,
    hydrogen_bonds_kernel, pi_stacking_kernel, salt_bridges_kernel, surface_contacts_kernel,
    water_bridges_kernel,
};
pub use kernel::{
    FrameKernelResult, MappedClosureStructureKernel, StructureKernel, mapped_structure_kernel,
    structure_kernel,
};
pub use native::{GovernedNativeError, native_contact_fraction_kernel};
pub use physical::{
    PhysicalKernelError, centre_of_mass_radial_distribution_kernel, coordination_numbers_kernel,
    density_map_kernel, leaflets_kernel, linear_density_kernel, pore_profile_kernel,
    radial_distribution_kernel,
};
pub use standalone::{
    StandaloneAnalysisError, governed_fragment_mapping, governed_polymer_statistics,
    governed_sugar_pucker,
};
pub use structure::{
    chain_interface_kernel, gnm_kernel, half_sphere_exposure_kernel, nucleic_torsions_kernel,
    secondary_structure_kernel,
};

#[cfg(test)]
#[path = "adapters_tests.rs"]
mod adapters_tests;

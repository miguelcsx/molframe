//! Contacts, interactions and secondary structure over a structure.
//!
//! Everything here reads positions and topology from a [`Structure`] and returns
//! deterministic, sorted findings. Neighbour searches go through the shared
//! spatial indices rather than an all-pairs scan, and any distance test is made
//! against a squared cutoff so the square root is only taken when a caller wants
//! a length reported.
//!
//! [`Structure`]: molframe_core::structure::Structure

#![forbid(unsafe_code)]

#[macro_use]
mod tables;

#[path = "anm.rs"]
mod anisotropic_network;
#[path = "pi_stacking.rs"]
mod aromatic_stacking;
#[path = "contacts.rs"]
mod atom_pairs;
#[path = "cation_pi.rs"]
mod cation_aromatic;
#[path = "interface.rs"]
mod chain_boundary;
#[path = "polymer.rs"]
mod chain_statistics;
#[path = "pore.rs"]
mod channel_profile;
mod chemistry;
#[path = "gnm.rs"]
mod elastic_network;
#[path = "dssp_binary.rs"]
mod external_secondary_structure;
#[path = "fragment.rs"]
mod fragment_mapping;
#[path = "hse.rs"]
mod half_sphere;
pub mod hbond;
#[path = "helical.rs"]
mod helix_geometry;
#[path = "leaflet.rs"]
mod membrane_layers;
#[path = "nmd.rs"]
mod mode_interchange;
#[path = "dynamics.rs"]
mod motion_statistics;
mod network;
#[path = "nucleic.rs"]
mod nucleic_torsion;
mod numeric;
#[path = "radial.rs"]
mod pair_distribution;
#[path = "base_pair.rs"]
mod paired_bases;
#[path = "governed/mod.rs"]
mod policy_execution;
#[path = "native.rs"]
mod reference_contacts;
#[path = "contact_map.rs"]
mod residue_contacts;
pub mod salt_bridge;
#[path = "dssp.rs"]
mod secondary_structure_assignment;
#[path = "density.rs"]
mod spatial_density;
mod stream_surface;
#[path = "pucker.rs"]
mod sugar_conformation;
pub mod surface_contacts;
#[path = "ensemble.rs"]
mod trajectory_ensemble;
mod vector_field;
#[path = "water_bridge.rs"]
mod water_mediation;

#[cfg(test)]
mod chemistry_test_support;

pub use anisotropic_network::{
    AnisotropicNetworkModel, AnmError, AnmOptions, anisotropic_network_model,
};
pub use aromatic_stacking::{
    PiStacking, PiStackingError, PiStackingOptions, PiStackingTable, StackingKind, pi_stacking,
};
pub use atom_pairs::{
    Contact, ContactTable, atom_contacts, atom_contacts_between,
    atom_contacts_between_with_spatial, visit_atom_contacts, visit_atom_contacts_between,
};
pub use cation_aromatic::{CationPi, CationPiError, CationPiOptions, CationPiTable, cation_pi};
pub use chain_boundary::{chain_interface, chain_interface_with_spatial};
pub use chain_statistics::{PolymerError, PolymerStatistics, polymer_statistics};
pub use channel_profile::{PoreError, PoreProfileOptions, PoreSample, pore_profile};
pub use elastic_network::{GaussianNetworkModel, GnmError, GnmOptions, gaussian_network_model};
pub use external_secondary_structure::{DsspBinaryError, DsspSegment, parse_dssp_output, run_dssp};
pub use fragment_mapping::{FragmentMappingError, FragmentMatch, FragmentReference, map_fragments};
pub use half_sphere::{HalfSphereExposure, HseError, half_sphere_exposure};
pub use hbond::{
    HydrogenBond, HydrogenBondError, HydrogenBondOptions, HydrogenBondTable, hydrogen_bonds,
};
pub use helix_geometry::{
    BaseFrame, HelicalError, HelicalOptions, HelicalParameters, helical_parameters, helical_steps,
};
pub use membrane_layers::{Leaflet, LeafletOptions, identify_leaflets};
pub use mode_interchange::{NmdError, NormalMode, NormalModeSet, read_nmd, write_nmd};
pub use motion_statistics::{
    DielectricOptions, DielectricResult, DynamicsError, WaterDynamicsOptions, WaterLag,
    dielectric_from_dipoles, water_dynamics,
};
pub use nucleic_torsion::{NucleicTorsionError, NucleicTorsions, nucleic_torsions};
pub use pair_distribution::{
    CentreGroup, CoordinationOptions, RadialBin, RadialDistributionOptions, RadialError,
    centre_of_mass_radial_distribution, coordination_numbers, radial_distribution,
};
pub use paired_bases::{BasePair, BasePairError, BasePairOptions, base_pairs};
pub use policy_execution::{
    AnalysisDescriptor, FrameKernelResult, GovernedAnalysisError, GovernedNativeError,
    GovernedStructureAnalysis, MappedClosureStructureKernel, PhysicalKernelError,
    StandaloneAnalysisError, StructureKernel, analyse_structure, analyse_trajectory,
    base_pairs_kernel, cation_pi_kernel, centre_of_mass_radial_distribution_kernel,
    chain_interface_kernel, contact_map_kernel, contacts_kernel, coordination_numbers_kernel,
    density_map_kernel, gnm_kernel, governed_dielectric_from_dipoles, governed_fragment_mapping,
    governed_helical_parameters, governed_helical_steps, governed_polymer_statistics,
    governed_sugar_pucker, governed_water_dynamics, half_sphere_exposure_kernel,
    hydrogen_bonds_kernel, leaflets_kernel, linear_density_kernel, mapped_structure_kernel,
    native_contact_fraction_kernel, nucleic_torsions_kernel, pi_stacking_kernel,
    pore_profile_kernel, radial_distribution_kernel, salt_bridges_kernel,
    secondary_structure_kernel, structure_kernel, surface_contacts_kernel, water_bridges_kernel,
};
pub use reference_contacts::{NativeContacts, NativeError, native_contact_fraction};
pub use residue_contacts::{ContactMap, ResidueContact, ResidueContactTable, residue_contact_map};
pub use salt_bridge::{SaltBridge, SaltBridgeTable, salt_bridges};
pub use secondary_structure_assignment::{
    DsspError, DsspOptions, SseKind, SseRecord, SseTable, secondary_structure,
};
pub use spatial_density::{
    CartesianAxis, DensityError, DensityGrid, DensityGridSpec, LinearDensityBin,
    LinearDensityOptions, density_map, linear_density,
};
pub use sugar_conformation::{Pucker, sugar_pucker};
pub use surface_contacts::{SurfaceContactOptions, surface_contacts};
pub use trajectory_ensemble::{
    Clustering, ConvergenceBlock, EnsembleDistanceMatrix, EnsembleGeometryError,
    EnsembleSimilarityError, EnsembleStatisticsError, FrameAlignment, GovernedEnsembleError,
    GroupVariance, HarmonicSimilarity, HarmonicSimilarityOptions, KMeans, KMeansError,
    KMeansOptions, Linkage, RemainderPolicy, agglomerative_clustering,
    analyse_agglomerative_clustering, analyse_block_convergence,
    analyse_cluster_population_similarity, analyse_dbscan_clustering,
    analyse_generalized_procrustes_mean, analyse_group_coordinate_variance,
    analyse_harmonic_ensemble_similarity, analyse_kmeans, analyse_medoid,
    analyse_pairwise_fitted_rmsd, analyse_pairwise_torus_distance, analyse_rmsd_to_reference,
    block_convergence, cluster_population_similarity, dbscan_clustering,
    generalized_procrustes_mean, group_coordinate_variance, harmonic_ensemble_similarity, kmeans,
    medoid, pairwise_fitted_rmsd, pairwise_torus_distance, rmsd_to_reference,
};
pub use vector_field::{
    StreamlineDirection, StreamlineOptions, VectorFieldError, VectorFieldGrid,
    integrate_streamlines,
};
pub use water_mediation::{WaterBridge, WaterBridgeOptions, WaterBridgeTable, water_bridges};

pub use stream_surface::{SasaStreamError, sasa_stream};

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

#[cfg(test)]
mod cif_fixture;
mod diagnostics;
mod numeric;
mod policy_execution;

mod elastic;
use elastic::{anisotropic_network, elastic_network, mode_interchange, network, vector_field};
mod interaction;
pub use interaction::hbond;
pub use interaction::salt_bridge;
pub use interaction::surface_contacts;
use interaction::{
    aromatic_rings, aromatic_stacking, atom_pairs, cation_aromatic, chain_boundary,
    channel_profile, chemistry, electrostatics, fragment_mapping, half_sphere, native_match,
    protonation, reference_contacts, residue_contacts, salt_bridge_groups, water_mediation,
};
mod nucleotide;
use nucleotide::{helix_geometry, nucleic_torsion, paired_bases, sugar_conformation};
mod trajectory;
use trajectory::{
    chain_statistics, membrane_layers, motion_statistics, pair_distribution, spatial_density,
    trajectory_ensemble,
};
mod secondary;
use secondary::{external_secondary_structure, secondary_structure_assignment};
mod area;
use area::{max_area_table, residue_sasa, sasa_atoms, stream_plan, stream_surface};

#[cfg(test)]
mod chemistry_test_support;
#[cfg(test)]
mod ring_test_support;

pub use anisotropic_network::{
    AnisotropicNetworkModel, AnmError, AnmOptions, anisotropic_network_model,
};
pub use aromatic_rings::{AromaticRing, MAXIMUM_RING_SIZE, aromatic_rings};
pub use aromatic_stacking::{
    PiStacking, PiStackingError, PiStackingGeometry, PiStackingOptions, PiStackingTable,
    StackingKind, pi_stacking, pi_stacking_with_geometry,
};
pub use atom_pairs::{
    Contact, ContactTable, atom_contacts, atom_contacts_between,
    atom_contacts_between_with_spatial, atom_contacts_periodic, visit_atom_contacts,
    visit_atom_contacts_between,
};
pub use cation_aromatic::{
    CationPi, CationPiError, CationPiOptions, CationPiTable, cation_pi, cation_pi_periodic,
};
pub use chain_boundary::{
    ChainNamespace, ChainSelector, InterfaceError, InterfaceOptions, chain_interface,
    chain_interface_with_options, chain_interface_with_spatial,
};
pub use chain_statistics::{PolymerError, PolymerStatistics, polymer_statistics};
pub use channel_profile::{PoreError, PoreProfileOptions, PoreSample, pore_profile};
pub use diagnostics::{governed_diagnostic, physical_diagnostic, standalone_diagnostic};
pub use elastic_network::{GaussianNetworkModel, GnmError, GnmOptions, gaussian_network_model};
pub use electrostatics::{
    CONTACT_POTENTIAL_CUTOFF, GridSpec, PotentialError, ScalarGrid, contact_potential,
    contact_potential_in,
};
pub use external_secondary_structure::{DsspBinaryError, DsspSegment, parse_dssp_output, run_dssp};
pub use fragment_mapping::{FragmentMappingError, FragmentMatch, FragmentReference, map_fragments};
pub use half_sphere::{HalfSphereExposure, HseError, half_sphere_exposure, protein_domain};
pub use hbond::{
    HeavyAtomHydrogenBond, HeavyAtomHydrogenBondOptions, HeavyAtomHydrogenBondTable,
    HeavyRoleTable, HydrogenBond, HydrogenBondCandidate, HydrogenBondCandidateOptions,
    HydrogenBondCandidateTable, HydrogenBondError, HydrogenBondOptions, HydrogenBondPolicy,
    HydrogenBondTable, heavy_atom_hydrogen_bonds, hydrogen_bond_candidates_over_segment,
    hydrogen_bonds, hydrogen_bonds_at_sample, hydrogen_bonds_with_policy,
};
pub use helix_geometry::{
    BaseFrame, HelicalError, HelicalOptions, HelicalParameters, MiddleFrameStep,
    helical_parameters, helical_steps, middle_frame_step, middle_frame_steps,
};
pub use max_area_table::MaxAreaTable;
pub use membrane_layers::{Leaflet, LeafletOptions, identify_leaflets};
pub use mode_interchange::{NmdError, NormalMode, NormalModeSet, read_nmd, write_nmd};
pub use molframe_surface::{SasaMethod, SasaOptions, SasaWordError};
pub use motion_statistics::{
    DielectricOptions, DielectricResult, DynamicsError, WaterDynamicsOptions, WaterLag,
    dielectric_from_dipoles, water_dynamics,
};
pub use native_match::match_atoms_by_name;
pub use nucleic_torsion::{
    NucleicTorsionError, NucleicTorsions, nucleic_torsions, nucleic_torsions_with_policy,
};
pub use pair_distribution::{
    CentreGroup, CoordinationOptions, RadialBin, RadialDistributionOptions, RadialError,
    centre_of_mass_radial_distribution, coordination_numbers, radial_distribution,
};
pub use paired_bases::{BasePair, BasePairError, BasePairOptions, WatsonCrickGeometry, base_pairs};
pub use policy_execution::{
    AnalysisDescriptor, DefinitionError, FrameKernelResult, FrameRecord, GovernedAnalysisError,
    GovernedNativeError, GovernedStructureAnalysis, MappedClosureStructureKernel,
    NativeCorrespondence, PhysicalKernelError, ResidueKernelError, StandaloneAnalysisError,
    StructureKernel, SurfaceSampling, analyse_models, analyse_structure, analyse_trajectory,
    base_pairs_kernel, cation_pi_kernel, centre_of_mass_radial_distribution_kernel,
    chain_interface_kernel, contact_map_kernel, contacts_kernel, coordination_numbers_kernel,
    definition_contacts_kernel, density_map_kernel, gnm_kernel, governed_dielectric_from_dipoles,
    governed_fragment_mapping, governed_helical_parameters, governed_helical_steps,
    governed_polymer_statistics, governed_sugar_pucker, governed_water_dynamics,
    half_sphere_exposure_kernel, hydrogen_bonds_kernel, leaflets_kernel, linear_density_kernel,
    mapped_structure_kernel, native_contact_fraction_kernel, native_contact_fraction_kernel_with,
    nucleic_torsions_kernel, pi_stacking_kernel, pore_profile_kernel, radial_distribution_kernel,
    residue_sasa_kernel, salt_bridges_kernel, sasa_kernel, secondary_structure_kernel,
    structure_kernel, surface_contacts_kernel, water_bridges_kernel,
};
pub use protonation::parent_residue;
pub use reference_contacts::{
    NativeContacts, NativeError, native_contact_fraction, native_contact_fraction_matched,
};
pub use residue_contacts::{
    ContactMap, ResidueContact, ResidueContactTable, residue_contact_map,
    residue_contact_map_periodic,
};
pub use residue_sasa::{ResidueArea, ResidueAreas, ResidueSasaError, residue_areas};
pub use salt_bridge::{
    SaltBridge, SaltBridgeError, SaltBridgeOptions, SaltBridgeTable, salt_bridges,
    salt_bridges_with_options,
};
pub use secondary_structure_assignment::{
    DsspError, DsspOptions, SseRecord, SseTable, secondary_structure,
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
pub use water_mediation::{
    WaterBridge, WaterBridgeConfidence, WaterBridgeMode, WaterBridgeOptions, WaterBridgePolicy,
    WaterBridgeReport, WaterBridgeTable, water_bridges, water_bridges_with_policy,
};

pub use stream_plan::{SasaSeries, SasaStreamPlan};
pub use stream_surface::{SasaStreamError, sasa_stream};

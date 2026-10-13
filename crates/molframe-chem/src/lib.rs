//! Deposited-component chemistry and versioned reference data.

#![forbid(unsafe_code)]

mod annotate;
mod aromaticity;
pub mod base_pair;
mod bonds;
mod carbohydrates;
mod coverage;
mod diagnostics;
mod element;
mod element_data;
mod equivalence;
mod grid;
mod ionic;
mod library;
mod metal;
mod model;
mod mol;
mod mol2;
mod molecule_structure;
mod numeric;
mod peoe;
mod provider;
mod radii_defaults;
mod radius_policy;
mod radius_policy_text;
mod radius_report;
mod radius_source;
mod radius_table;
mod roles;
mod secondary;
mod side_chain;
mod smarts;
mod smarts_match;
mod smarts_parse;
mod standard_bonds;
mod standard_components;
mod standard_roles;
mod table_identity;

pub use annotate::{
    ChemistryProvenance, ChemistryReport, FormalChargeSource, HBOND_WEAK_ACCEPTOR_ANNOTATION,
    PolymerLinkPolicy, PolymerLinkRule, apply_component_chemistry,
};
pub use bonds::{covalent_pair, perceive_bonds, perceive_bonds_in};
pub use carbohydrates::{
    CarbohydrateLink, CarbohydrateOptions, CarbohydrateReport, Monosaccharide, RingGeometry,
    SnfgShape, SnfgSymbol, carbohydrates, snfg_symbol,
};
pub use coverage::{ComponentCoverage, component_coverage};
pub use element::{
    ElementProperties, RadiusSet, RadiusTable, atom_masses, atom_radii, element_properties,
    vdw_radius,
};
pub use equivalence::{
    AutomorphismLimit, EquivalenceCache, EquivalenceClasses, automorphisms, equivalence_classes,
};
pub use ionic::{IonicRadius, IonicSpin, ionic_radii};
pub use library::{
    FINGERPRINT_BITS, Fingerprint, LigandDescriptors, RuleOfFive, cluster, describe, fingerprint,
    tanimoto,
};
pub use metal::{METAL_SITES, MetalRadiusKind, MetalSite, metal_radius, metal_site};
pub use model::{
    BondStereo, Component, ComponentAtom, ComponentBond, ComponentKind, PolymerAtomRole,
    StereoConfiguration,
};
pub use mol::{
    MolAtom, MolAtomMetadata, MolBond, MolBondMetadata, MolError, MolRecord, MolVersion, Molecule,
    SdfProperty, parse_mol_record, parse_sdf_records, write_mol, write_sdf,
};
pub use mol2::{
    Mol2AtomMetadata, Mol2BondMetadata, Mol2Error, Mol2Record, Mol2Section, parse_mol2_record,
    write_mol2,
};
pub use molecule_structure::{
    mol_record_to_structure, mol2_record_to_structure, molecule_to_structure,
    small_cif_to_structure, structure_to_molecule,
};
pub use peoe::{
    ChargeSource, PartialChargeError, PartialCharges, PeoeAtom, PeoeAtomType, PeoeBond, PeoeError,
    PeoeOptions, PeoeParameterProfile, component_peoe_charges, partial_charges, peoe_charges,
};
pub use provider::{
    CORE_DICTIONARY_VERSION, CifProvider, ComponentProvider, CoreProvider,
    EMBEDDED_COMPONENTS_EXTENSION, EmbeddedProvider, LayeredProvider, LazyCifProvider,
    MemoryProvider, embedded_components, read_ccd, resolution_stack, take_embedded_components,
};
pub use radii_defaults::{
    ALVAREZ as ALVAREZ_RADII, AMBER_UNITED as AMBER_UNITED_RADII, BONDI as BONDI_RADII,
    CHARMM as CHARMM_RADII, METAL_SITE_TABLE, PYYKKO_COVALENT as PYYKKO_COVALENT_RADII,
    SHANNON_IONIC as SHANNON_IONIC_RADII,
};
pub use radius_policy::{
    HeteroPolicy, MetalFallback, RadiusPolicy, ResidueClass, classify_residue, included_atoms,
    policy_radii, radius_with_fallback, source_radius_with_fallback,
};
pub use radius_policy_text::PolicyWordError;
pub use radius_report::AppliedRadiusPolicy;
pub use radius_source::RadiusSource;
pub use radius_table::{CustomRadii, RadiusTableError};
pub use roles::{
    PolymerRoleProfile, PolymerRoleReport, PolymerRoleRule, apply_polymer_role_profile,
};
pub use secondary::{
    DsspBackbone, DsspOptions, InvalidDsspOptions, assign_secondary_structure,
    assign_secondary_structure_in, dssp_from_backbones,
};
pub use side_chain::{SideChainDefinition, SideChainRoles, side_chain_definition};
pub use smarts::{SmartsDataError, SmartsError, SmartsMatch, SmartsPattern};
pub use standard_bonds::{
    is_amino_acid_component, is_nucleotide_component, is_water_component, standard_bond_order,
};
pub use standard_components::{annotate_standard_components, standard_component_kind};
pub use standard_roles::{
    STANDARD_PROTEIN_ROLE_PROFILE, annotate_standard_protein_roles, standard_protein_atom_role,
};
pub use table_identity::TableIdentity;

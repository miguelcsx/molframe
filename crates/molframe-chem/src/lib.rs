//! Deposited-component chemistry and versioned reference data.
//!
//! Reference tables are selected through a [`TableStore`], which consults three
//! layers in turn: tables and a profile registered in memory
//! ([`TableStore::insert_profile`], [`TableStore::insert_table`]), then the
//! profile and cache directory on disk, then the tables the `embedded-tables`
//! feature compiles in (the bonding reach alone). Without that feature nothing
//! is compiled in.

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
mod formats;
mod grid;
mod library;
mod model;
mod numeric;
mod pattern;
mod peoe;
mod provider;
mod radii;
mod roles;
mod secondary;
mod side_chain;
mod standard;
mod tables;
#[cfg(test)]
mod test_data;

use formats::{mol, mol2, molecule_structure};
use pattern::{smarts, smarts_match, smarts_parse};
use radii::{
    bond_reach, element_radii, ionic, metal, radius_policy, radius_policy_text, radius_report,
    radius_selection, radius_source,
};
use standard::{standard_bonds, standard_components, standard_roles};
use tables::{
    table_bundled, table_bytes, table_identity, table_memory, table_provider, table_request,
    table_row, table_store, table_text,
};

pub use annotate::{
    ChemistryProvenance, ChemistryReport, FormalChargeSource, HBOND_WEAK_ACCEPTOR_ANNOTATION,
    PolymerLinkPolicy, PolymerLinkRule, apply_component_chemistry,
};
pub use bond_reach::BondReach;
pub use bonds::{covalent_pair, perceive_bonds, perceive_bonds_in};
pub use carbohydrates::{
    CarbohydrateLink, CarbohydrateOptions, CarbohydrateReport, Monosaccharide, RingGeometry,
    SnfgShape, SnfgSymbol, carbohydrates, snfg_symbol,
};
pub use coverage::{ComponentCoverage, component_coverage};
pub use element::{ElementProperties, atom_masses, atom_radii, element_properties};
pub use element_radii::{Covalent, CovalentRadii, ElementRadii, RadiusRole, Vdw, VdwRadii};
pub use equivalence::{
    AutomorphismLimit, EquivalenceCache, EquivalenceClasses, automorphisms, equivalence_classes,
};
pub use ionic::{IonicRadii, IonicRadius, IonicSpin};
pub use library::{
    FINGERPRINT_BITS, Fingerprint, LigandDescriptors, RuleOfFive, cluster, describe, fingerprint,
    tanimoto,
};
pub use metal::{MetalRadiusKind, MetalSite, MetalSites, MetalTables};
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
    CifProvider, ComponentDictionary, ComponentProvider, EMBEDDED_COMPONENTS_EXTENSION,
    EmbeddedProvider, LayeredProvider, LazyCifProvider, MemoryProvider, embedded_components,
    profile_dictionary, read_ccd, resolution_stack, take_embedded_components,
};
pub use radius_policy::{
    HeteroPolicy, MetalFallback, RadiusPolicy, ResidueClass, classify_residue, included_atoms,
    policy_radii,
};
pub use radius_policy_text::PolicyWordError;
pub use radius_report::AppliedRadiusPolicy;
pub use radius_selection::RadiusSelection;
pub use radius_source::RadiusSource;
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
pub use table_bytes::TableBytes;
pub use table_identity::TableIdentity;
pub use table_provider::{
    Cache, InMemory, MemoryFile, Provided, Source, TableError, TableKind, TableProvider, UserFile,
};
pub use table_request::TableRequest;
pub use table_row::Row;
pub use table_store::{Selected, TableStore};
pub use table_text::{
    IdentityOnly, TABLE_FORMAT_VERSION, TableLines, TableTextError, join_binary, read_table,
    split_binary, text_of,
};

//! Deposited-component chemistry and versioned reference data.

#![forbid(unsafe_code)]

mod annotate;
mod coverage;
mod element;
mod element_data;
mod equivalence;
mod ionic;
mod model;
mod mol;
mod mol2;
mod numeric;
mod peoe;
mod provider;
mod roles;
mod side_chain;
mod smarts;
mod smarts_match;
mod smarts_parse;

pub use annotate::{
    ChemistryReport, PolymerLinkPolicy, PolymerLinkRule, apply_component_chemistry,
};
pub use coverage::{ComponentCoverage, component_coverage};
pub use element::{ElementProperties, RadiusSet, RadiusTable, element_properties, vdw_radius};
pub use equivalence::{
    AutomorphismLimit, EquivalenceCache, EquivalenceClasses, automorphisms, equivalence_classes,
};
pub use ionic::{IonicRadius, IonicSpin, ionic_radii};
pub use model::{
    Component, ComponentAtom, ComponentBond, ComponentKind, PolymerAtomRole, StereoConfiguration,
};
pub use mol::{
    MolAtom, MolAtomMetadata, MolBond, MolBondMetadata, MolError, MolRecord, MolVersion, Molecule,
    SdfProperty, parse_mol_record, parse_sdf_records, write_mol, write_sdf,
};
pub use mol2::{
    Mol2AtomMetadata, Mol2BondMetadata, Mol2Error, Mol2Record, Mol2Section, parse_mol2_record,
    write_mol2,
};
pub use peoe::{
    PeoeAtom, PeoeAtomType, PeoeBond, PeoeError, PeoeOptions, PeoeParameterProfile,
    component_peoe_charges, peoe_charges,
};
pub use provider::{CifProvider, ComponentProvider, MemoryProvider, read_ccd};
pub use roles::{
    PolymerRoleProfile, PolymerRoleReport, PolymerRoleRule, apply_polymer_role_profile,
};
pub use side_chain::{SideChainDefinition, SideChainRoles, side_chain_definition};
pub use smarts::{SmartsDataError, SmartsError, SmartsMatch, SmartsPattern};

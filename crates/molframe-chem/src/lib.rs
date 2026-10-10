//! Deposited-component chemistry and versioned reference data.

#![forbid(unsafe_code)]

mod annotate;
mod aromaticity;
mod bonds;
mod carbohydrates;
mod coverage;
mod diagnostics;
mod element;
mod element_data;
mod equivalence;
mod grid;
mod ionic;
mod model;
mod mol;
mod mol2;
mod molecule_structure;
mod numeric;
mod peoe;
mod provider;
mod roles;
mod secondary;
mod side_chain;
mod smarts;
mod smarts_match;
mod smarts_parse;
mod standard_bonds;
mod standard_components;

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
pub use molecule_structure::{
    mol_record_to_structure, mol2_record_to_structure, molecule_to_structure,
    small_cif_to_structure, structure_to_molecule,
};
pub use peoe::{
    ChargeSource, PartialChargeError, PartialCharges, PeoeAtom, PeoeAtomType, PeoeBond, PeoeError,
    PeoeOptions, PeoeParameterProfile, component_peoe_charges, partial_charges, peoe_charges,
};
pub use provider::{CifProvider, ComponentProvider, MemoryProvider, read_ccd};
pub use roles::{
    PolymerRoleProfile, PolymerRoleReport, PolymerRoleRule, apply_polymer_role_profile,
};
pub use secondary::{
    DsspBackbone, DsspOptions, InvalidDsspOptions, assign_secondary_structure, dssp_from_backbones,
};
pub use side_chain::{SideChainDefinition, SideChainRoles, side_chain_definition};
pub use smarts::{SmartsDataError, SmartsError, SmartsMatch, SmartsPattern};
pub use standard_bonds::{
    is_amino_acid_component, is_nucleotide_component, is_water_component, standard_bond_order,
};
pub use standard_components::{annotate_standard_components, standard_component_kind};

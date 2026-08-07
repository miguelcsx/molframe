//! Deposited-component chemistry and versioned reference data.

#![forbid(unsafe_code)]

mod annotate;
mod coverage;
mod element;
mod element_data;
mod equivalence;
mod ionic;
mod model;
mod provider;
mod side_chain;

pub use annotate::{ChemistryReport, apply_component_chemistry};
pub use coverage::{ComponentCoverage, component_coverage};
pub use element::{ElementProperties, RadiusSet, RadiusTable, element_properties, vdw_radius};
pub use equivalence::{EquivalenceCache, EquivalenceClasses, equivalence_classes};
pub use ionic::{IonicRadius, IonicSpin, ionic_radii};
pub use model::{Component, ComponentAtom, ComponentBond, ComponentKind, StereoConfiguration};
pub use provider::{CifProvider, ComponentProvider, MemoryProvider, read_ccd};
pub use side_chain::{SideChainDefinition, side_chain_definition};

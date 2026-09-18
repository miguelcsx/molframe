//! Hybridisation-aware partial equalisation of orbital electronegativity.
//!
//! Component topology is perceived once, before iteration. Charge transfer is
//! then `O(passes * bonds)` and allocates no memory inside the bond loop.

mod calculate;
mod parameters;
mod perceive;
mod types;

pub use calculate::{component_peoe_charges, peoe_charges};
pub use types::{PeoeAtom, PeoeAtomType, PeoeBond, PeoeError, PeoeOptions, PeoeParameterProfile};

#[cfg(test)]
#[path = "peoe_tests.rs"]
mod tests;

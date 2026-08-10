//! Interactive Molecular Dynamics streaming protocol.

mod client;
mod error;
mod model;
mod wire;

pub use client::ImdClient;
pub use error::ImdError;
pub use model::{
    ImdConnectionOptions, ImdEnergies, ImdForce, ImdLimits, ImdMessage, ImdPeerEndian,
};

#[cfg(test)]
#[path = "imd_tests.rs"]
mod tests;

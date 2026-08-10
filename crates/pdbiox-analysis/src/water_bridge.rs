//! Water-mediated networks assembled from fully perceived hydrogen bonds.

use std::collections::BTreeMap;

use pdbiox_core::index::AtomIndex;
use pdbiox_core::structure::Structure;

use crate::hbond::{HydrogenBondError, HydrogenBondOptions, hydrogen_bonds};

/// Explicit water-bridge detection policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterBridgeOptions {
    /// Geometry and periodic policy for each constituent hydrogen bond.
    pub hydrogen_bonds: HydrogenBondOptions,
}

/// A solvent heavy atom bridging two non-solvent partners.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WaterBridge {
    /// CCD-classified solvent donor or acceptor atom.
    pub water: AtomIndex,
    /// Lower-indexed non-solvent partner.
    pub first: AtomIndex,
    /// Higher-indexed non-solvent partner.
    pub second: AtomIndex,
}

/// Builds water bridges from two or more oriented hydrogen bonds to one solvent.
///
/// Component identity and donor/acceptor roles come exclusively from CCD
/// annotations. Heavy-atom and residue-name fallbacks are not performed.
///
/// # Errors
///
/// Propagates hydrogen-bond chemistry, geometry, periodic and spatial errors.
pub fn water_bridges(
    structure: &Structure,
    options: WaterBridgeOptions,
) -> Result<Vec<WaterBridge>, HydrogenBondError> {
    let bonds = hydrogen_bonds(structure, options.hydrogen_bonds)?;
    let mut partners: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for bond in bonds {
        let donor_water = crate::chemistry::component_kind(structure, bond.donor.get())
            == Some(pdbiox_chem::ComponentKind::Solvent);
        let acceptor_water = crate::chemistry::component_kind(structure, bond.acceptor.get())
            == Some(pdbiox_chem::ComponentKind::Solvent);
        match (donor_water, acceptor_water) {
            (true, false) => partners
                .entry(bond.donor.get())
                .or_default()
                .push(bond.acceptor.get()),
            (false, true) => partners
                .entry(bond.acceptor.get())
                .or_default()
                .push(bond.donor.get()),
            (true, true) | (false, false) => {}
        }
    }
    let mut bridges = Vec::new();
    for (water, mut atoms) in partners {
        atoms.sort_unstable();
        atoms.dedup();
        for first in 0..atoms.len() {
            for second in first + 1..atoms.len() {
                bridges.push(WaterBridge {
                    water: AtomIndex::new(water),
                    first: AtomIndex::new(atoms[first]),
                    second: AtomIndex::new(atoms[second]),
                });
            }
        }
    }
    Ok(bridges)
}

#[cfg(test)]
#[path = "water_bridge_tests.rs"]
mod tests;

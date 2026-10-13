//! One declared radius convention for the metal ions of structural biology.
//!
//! Bond perception, solvent-accessible surface and steric clash detection all
//! meet the same ions and need three different radii for them: a reach inside
//! which a ligand counts as coordinated, a van der Waals radius that sizes the
//! ion as a solvent obstacle, and an ionic radius that sizes a coordinated
//! cation as a steric contact. Keeping the three in one table, keyed by
//! element, is what makes them agree with each other.
//!
//! Each metal carries the oxidation state and coordination it most often has in
//! deposited macromolecules; the ionic radius is Shannon's effective radius for
//! exactly that pair, read from the bundled Shannon catalogue rather than
//! copied here, so a change of catalogue moves every consumer together. The van
//! der Waals radius is Alvarez's (2013), because Bondi's set has no value for
//! most of these elements and the library never fills one set from another
//! silently; the fallback is a named, explicit choice.
//!
//! Lookups are constant time on a short static table and allocate nothing.

use crate::ionic::{IonicSpin, ionic_radii};
use crate::{RadiusSet, vdw_radius};
use molframe_core::Element;

/// Which radius of a metal ion a consumer wants.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MetalRadiusKind {
    /// Alvarez van der Waals radius, the size of the ion as a solvent obstacle.
    Vdw,
    /// Shannon effective ionic radius at the metal's typical charge and
    /// coordination, the size of a coordinated cation as a steric contact.
    Ionic,
    /// Single-bond covalent radius.
    Covalent,
}

/// The oxidation state and coordination a metal usually has in macromolecules,
/// and the reach inside which a ligand is perceived as coordinated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MetalSite {
    /// Element symbol.
    pub symbol: &'static str,
    /// Typical formal charge.
    pub charge: i8,
    /// Typical coordination in Shannon's notation.
    pub coordination: &'static str,
    /// Spin state the typical site takes.
    pub spin: IonicSpin,
    /// Bonding reach in ångström: a ligand with reach `r` is coordinated up to
    /// `(reach + r) / 1.95`.
    pub reach: f32,
}

const fn site(
    symbol: &'static str,
    charge: i8,
    coordination: &'static str,
    spin: IonicSpin,
    reach: f32,
) -> MetalSite {
    MetalSite {
        symbol,
        charge,
        coordination,
        spin,
        reach,
    }
}

const NONE: IonicSpin = IonicSpin::Unspecified;
const HIGH: IonicSpin = IonicSpin::High;

/// Every metal the library perceives coordination for, in atomic-number order.
///
/// A heme iron sits 1.9 to 2.3 Å from its pyrrole nitrogens and a calcium 2.3
/// to 2.6 Å from its oxygens, so one shared reach cannot serve; each value is
/// chosen so the longest ordinary coordination distance to an N, O or S donor
/// still falls inside the window while a non-coordinated neighbour does not.
pub const METAL_SITES: [MetalSite; 20] = [
    site("Li", 1, "VI", NONE, 2.80),
    site("Na", 1, "VI", NONE, 3.30),
    site("Mg", 2, "VI", NONE, 3.00),
    site("Al", 3, "VI", NONE, 2.80),
    site("K", 1, "VI", NONE, 3.90),
    site("Ca", 2, "VI", NONE, 3.40),
    site("V", 3, "VI", NONE, 2.90),
    site("Cr", 3, "VI", NONE, 2.90),
    site("Mn", 2, "VI", HIGH, 2.90),
    site("Fe", 2, "VI", HIGH, 2.90),
    site("Co", 2, "VI", HIGH, 2.90),
    site("Ni", 2, "VI", NONE, 2.90),
    site("Cu", 2, "V", NONE, 2.90),
    site("Zn", 2, "IV", NONE, 2.85),
    site("Mo", 6, "VI", NONE, 3.10),
    site("Ag", 1, "II", NONE, 3.10),
    site("Cd", 2, "VI", NONE, 3.10),
    site("W", 6, "VI", NONE, 3.10),
    site("Au", 3, "IVSQ", NONE, 2.90),
    site("Hg", 2, "II", NONE, 3.00),
];

/// The typical site of a metal, or `None` for any other element.
#[must_use]
pub fn metal_site(symbol: &str) -> Option<&'static MetalSite> {
    METAL_SITES.iter().find(|site| site.symbol == symbol)
}

pub(crate) fn is_metal(symbol: &str) -> bool {
    metal_site(symbol).is_some()
}

/// A radius of a metal ion in ångström, or `None` for an element that is not a
/// perceived metal or has no value of that kind.
#[must_use]
pub fn metal_radius(element: Element, kind: MetalRadiusKind) -> Option<f32> {
    let site = metal_site(element.symbol())?;
    match kind {
        MetalRadiusKind::Vdw => vdw_radius(element, RadiusSet::Alvarez),
        MetalRadiusKind::Covalent => crate::element_properties(element)?.covalent_radius,
        MetalRadiusKind::Ionic => ionic_radii(element)
            .ok()?
            .iter()
            .find(|entry| {
                entry.charge == site.charge
                    && entry.coordination.as_ref() == site.coordination
                    && entry.spin == site.spin
            })
            .map(|entry| entry.ionic_radius),
    }
}

#[cfg(test)]
#[path = "metal_tests.rs"]
mod tests;

//! Definitional properties of chemical elements: weight, period, group, valence
//! electrons and electronegativity.
//!
//! These are the periodic table itself, fixed by definition (IUPAC weights, the
//! ground-state configuration, the Pauling scale), not a scientific choice a
//! user would swap, so they are the one place data is compiled in. Everything
//! that is a choice (van der Waals, covalent and ionic radii, bonding reach, metal
//! sites) is a table selected by name. Lookup is constant time by atomic number
//! and allocates nothing.

use crate::ElementRadiusSource;
use crate::element_data::ELEMENTS;
use molframe_core::Element;

/// Scalar reference properties of one element.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElementProperties {
    /// Standard atomic weight, or the conventional mass number for elements
    /// without a standard atomic weight.
    pub atomic_weight: f64,
    /// Pauling electronegativity where defined.
    pub electronegativity: Option<f32>,
    /// Valence-electron count from the ground-state electronic configuration.
    pub valence_electrons: u8,
    /// Period in the periodic table.
    pub period: u8,
    /// Group in the periodic table, absent for the f-block convention.
    pub group: Option<u8>,
}

#[derive(Clone, Copy)]
pub(super) struct RawElement {
    electronegativity: Option<f32>,
    period: u8,
    group: u8,
    valence_electrons: u8,
}

impl RawElement {
    pub(super) const EMPTY: Self = Self::new(None, 0, 0, 0);

    pub(super) const fn new(
        electronegativity: Option<f32>,
        period: u8,
        group: u8,
        valence_electrons: u8,
    ) -> Self {
        Self {
            electronegativity,
            period,
            group,
            valence_electrons,
        }
    }
}

/// Returns reference properties by atomic number.
#[must_use]
pub fn element_properties(element: Element) -> Option<ElementProperties> {
    let z = element.atomic_number() as usize;
    if z == 0 {
        return None;
    }
    let raw = *ELEMENTS.get(z)?;
    Some(ElementProperties {
        atomic_weight: *ATOMIC_WEIGHTS.get(z)?,
        electronegativity: raw.electronegativity,
        valence_electrons: raw.valence_electrons,
        period: raw.period,
        group: (raw.group != 0).then_some(raw.group),
    })
}

/// The radius of every atom of `structure`, indexed by atom, from `source`; `NaN`
/// marks an atom whose element is unknown or has no value in the table.
///
/// Costs one table lookup per atom and one allocation for the result.
#[must_use]
pub fn atom_radii<R>(
    structure: &molframe_core::structure::Structure,
    source: &ElementRadiusSource<R>,
) -> Vec<f32> {
    let mut radii = vec![f32::NAN; structure.atom_count() as usize];
    for atom in structure.data().atoms() {
        let radius = atom.element().and_then(|element| source.radius(element));
        if let (Some(radius), Some(slot)) = (radius, radii.get_mut(atom.index().as_usize())) {
            *slot = radius;
        }
    }
    radii
}

/// The atomic weight of every atom of `structure`, indexed by atom, in daltons;
/// `NaN` marks an atom whose element is unknown.
#[must_use]
pub fn atom_masses(structure: &molframe_core::structure::Structure) -> Vec<f64> {
    let mut masses = vec![f64::NAN; structure.atom_count() as usize];
    for atom in structure.data().atoms() {
        let weight = atom
            .element()
            .and_then(element_properties)
            .map(|properties| properties.atomic_weight);
        if let (Some(weight), Some(slot)) = (weight, masses.get_mut(atom.index().as_usize())) {
            *slot = weight;
        }
    }
    masses
}

const ATOMIC_WEIGHTS: [f64; 119] = [
    0.0,
    1.008,
    4.002_602,
    6.94,
    9.012_183_1,
    10.81,
    12.011,
    14.007,
    15.999,
    18.998_403_163,
    20.1797,
    22.989_769_28,
    24.305,
    26.981_538_5,
    28.085,
    30.973_761_998,
    32.06,
    35.45,
    39.948,
    39.0983,
    40.078,
    44.955_908,
    47.867,
    50.9415,
    51.9961,
    54.938_044,
    55.845,
    58.933_194,
    58.6934,
    63.546,
    65.38,
    69.723,
    72.630,
    74.921_595,
    78.971,
    79.904,
    83.798,
    85.4678,
    87.62,
    88.905_84,
    91.224,
    92.906_37,
    95.95,
    98.0,
    101.07,
    102.905_50,
    106.42,
    107.8682,
    112.414,
    114.818,
    118.710,
    121.760,
    127.60,
    126.904_47,
    131.293,
    132.905_451_96,
    137.327,
    138.905_47,
    140.116,
    140.907_66,
    144.242,
    145.0,
    150.36,
    151.964,
    157.25,
    158.925_35,
    162.500,
    164.930_33,
    167.259,
    168.934_22,
    173.045,
    174.9668,
    178.49,
    180.947_88,
    183.84,
    186.207,
    190.23,
    192.217,
    195.084,
    196.966_569,
    200.592,
    204.38,
    207.2,
    208.980_40,
    209.0,
    210.0,
    222.0,
    223.0,
    226.0,
    227.0,
    232.0377,
    231.035_88,
    238.028_91,
    237.0,
    244.0,
    243.0,
    247.0,
    247.0,
    251.0,
    252.0,
    257.0,
    258.0,
    259.0,
    266.0,
    267.0,
    268.0,
    269.0,
    270.0,
    277.0,
    278.0,
    281.0,
    282.0,
    285.0,
    286.0,
    289.0,
    290.0,
    293.0,
    294.0,
    294.0,
];

#[cfg(test)]
#[path = "element_tests.rs"]
mod tests;

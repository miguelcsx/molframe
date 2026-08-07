//! Versioned properties of chemical elements.
//!
//! Lookup is constant time by atomic number and allocates nothing. Missing
//! values remain absent; the library never substitutes one radius convention
//! for another because that would change geometric results silently.

use crate::element_data::ELEMENTS;
use pdbiox_core::Element;

/// Scalar reference properties of one element.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElementProperties {
    /// Standard atomic weight, or the conventional mass number for elements
    /// without a standard atomic weight.
    pub atomic_weight: f64,
    /// Single-bond covalent radius in ångström.
    pub covalent_radius: Option<f32>,
    /// Pauling electronegativity where defined.
    pub electronegativity: Option<f32>,
    /// Valence-electron count from the ground-state electronic configuration.
    pub valence_electrons: u8,
    /// Period in the periodic table.
    pub period: u8,
    /// Group in the periodic table, absent for the f-block convention.
    pub group: Option<u8>,
}

/// A named van der Waals radius convention.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RadiusSet {
    /// Bondi's widely used crystallographic values.
    Bondi,
    /// AMBER united-atom values.
    AmberUnited,
    /// CHARMM force-field values.
    Charmm,
    /// Alvarez's crystallographic compilation.
    Alvarez,
}

/// Identity and version of a radius table.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RadiusTable {
    /// Stable machine-readable name.
    pub name: &'static str,
    /// Source/version label recorded in provenance.
    pub version: &'static str,
}

#[derive(Clone, Copy)]
pub(super) struct RawElement {
    covalent_radius: f32,
    electronegativity: Option<f32>,
    period: u8,
    group: u8,
    valence_electrons: u8,
    bondi: Option<f32>,
    alvarez: Option<f32>,
}

impl RawElement {
    pub(super) const EMPTY: Self = Self::new(0.0, None, 0, 0, 0, None, None);

    pub(super) const fn new(
        covalent_radius: f32,
        electronegativity: Option<f32>,
        period: u8,
        group: u8,
        valence_electrons: u8,
        bondi: Option<f32>,
        alvarez: Option<f32>,
    ) -> Self {
        Self {
            covalent_radius,
            electronegativity,
            period,
            group,
            valence_electrons,
            bondi,
            alvarez,
        }
    }
}

impl RadiusSet {
    /// The exact table identity represented by this value.
    #[must_use]
    pub const fn table(self) -> RadiusTable {
        match self {
            Self::Bondi => RadiusTable {
                name: "bondi",
                version: "1964",
            },
            Self::AmberUnited => RadiusTable {
                name: "amber_united",
                version: "weiner-1984",
            },
            Self::Charmm => RadiusTable {
                name: "charmm",
                version: "c36m-2021",
            },
            Self::Alvarez => RadiusTable {
                name: "alvarez",
                version: "2013",
            },
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
        covalent_radius: Some(raw.covalent_radius),
        electronegativity: raw.electronegativity,
        valence_electrons: raw.valence_electrons,
        period: raw.period,
        group: (raw.group != 0).then_some(raw.group),
    })
}

/// Returns a van der Waals radius in ångström from exactly the requested set.
#[must_use]
pub fn vdw_radius(element: Element, set: RadiusSet) -> Option<f32> {
    let z = element.atomic_number() as usize;
    match set {
        RadiusSet::Bondi => ELEMENTS.get(z).and_then(|element| element.bondi),
        RadiusSet::AmberUnited => amber(element),
        RadiusSet::Charmm => charmm(element),
        RadiusSet::Alvarez => ELEMENTS.get(z).and_then(|element| element.alvarez),
    }
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

fn amber(element: Element) -> Option<f32> {
    Some(match element.symbol() {
        "H" => 1.0,
        "C" | "N" => 1.8,
        "O" => 1.6,
        "F" => 1.55,
        "P" => 2.1,
        "S" | "Cl" => 2.0,
        "W" => 1.4,
        _ => return None,
    })
}

fn charmm(element: Element) -> Option<f32> {
    Some(match element.symbol() {
        "H" => 1.25,
        "C" | "S" => 2.0,
        "N" => 1.85,
        "O" => 1.70,
        "F" => 1.47,
        "P" => 2.15,
        "Cl" => 1.75,
        _ => return None,
    })
}

#[cfg(test)]
#[path = "element_tests.rs"]
mod tests;

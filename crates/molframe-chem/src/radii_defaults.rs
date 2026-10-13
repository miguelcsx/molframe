//! Default radius tables and the identity each one declares.
//!
//! This is the one place that says where the library's own radii come from.
//! The Bondi, Alvarez and covalent columns live in the committed element table
//! (`element_data.rs`, derived from mendeleev 1.1.0); the two force-field sets
//! are written out here. The algorithms that consume radii never name a value:
//! they receive a [`RadiusSource`](crate::RadiusSource), of which these tables
//! are only the defaults, so a caller's own table replaces them without any
//! change to the algorithm.

use crate::TableIdentity;

/// Bondi van der Waals radii.
pub const BONDI: TableIdentity = TableIdentity::builtin(
    "bondi",
    "1964",
    "Bondi A (1964) J. Phys. Chem. 68, 441, as distributed in mendeleev 1.1.0",
    "MIT (mendeleev 1.1.0 distribution)",
);

/// Alvarez van der Waals radii.
pub const ALVAREZ: TableIdentity = TableIdentity::builtin(
    "alvarez",
    "2013",
    "Alvarez S (2013) Dalton Trans. 42, 8617, as distributed in mendeleev 1.1.0",
    "MIT (mendeleev 1.1.0 distribution)",
);

/// AMBER united-atom radii.
pub const AMBER_UNITED: TableIdentity = TableIdentity::builtin(
    "amber_united",
    "weiner-1984",
    "Weiner SJ et al. (1984) J. Am. Chem. Soc. 106, 765, transcribed in this crate",
    "MIT (this repository)",
);

/// CHARMM radii.
pub const CHARMM: TableIdentity = TableIdentity::builtin(
    "charmm",
    "c36m-2021",
    "CHARMM36m van der Waals radii, transcribed in this crate",
    "MIT (this repository)",
);

/// Single-bond covalent radii.
pub const PYYKKO_COVALENT: TableIdentity = TableIdentity::builtin(
    "pyykko_covalent",
    "2009",
    "Pyykko P, Atsumi M (2009) Chem. Eur. J. 15, 186, as distributed in mendeleev 1.1.0",
    "MIT (mendeleev 1.1.0 distribution)",
);

/// Shannon effective ionic radii.
pub const SHANNON_IONIC: TableIdentity = TableIdentity::builtin(
    "shannon_ionic",
    "1976",
    "Shannon RD (1976) Acta Cryst. A32, 751, bundled as data/ionic-radii.bin",
    "MIT (this repository)",
);

/// The typical charge, coordination and bonding reach of each perceived metal.
pub const METAL_SITE_TABLE: TableIdentity = TableIdentity::builtin(
    "metal_sites",
    "1",
    "typical macromolecular oxidation state and coordination per metal, curated in this crate",
    "MIT (this repository)",
);

/// AMBER united-atom radii by element symbol, in ångström.
pub(crate) const AMBER_UNITED_RADII: &[(&str, f32)] = &[
    ("H", 1.0),
    ("C", 1.8),
    ("N", 1.8),
    ("O", 1.6),
    ("F", 1.55),
    ("P", 2.1),
    ("S", 2.0),
    ("Cl", 2.0),
    ("W", 1.4),
];

/// CHARMM radii by element symbol, in ångström.
pub(crate) const CHARMM_RADII: &[(&str, f32)] = &[
    ("H", 1.25),
    ("C", 2.0),
    ("S", 2.0),
    ("N", 1.85),
    ("O", 1.70),
    ("F", 1.47),
    ("P", 2.15),
    ("Cl", 1.75),
];

pub(crate) fn lookup(table: &[(&str, f32)], symbol: &str) -> Option<f32> {
    table
        .iter()
        .find(|(name, _)| *name == symbol)
        .map(|(_, radius)| *radius)
}

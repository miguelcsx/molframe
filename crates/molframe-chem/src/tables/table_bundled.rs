//! The tables compiled into the library, when the `embedded-tables` feature asks.
//!
//! A clean machine, or a target without a filesystem (a browser), has no
//! profile to select a table from, so without help the common flows refuse:
//! a read perceives no bonds, and bond inference, surface area and clash
//! scoring raise a missing-table error. With the feature, the store's last
//! layer is the minimum those default flows read, and a profile naming it:
//!
//! - `bond-reach: macromolecular-reach`, the reach read-time perception and
//!   crystal links use;
//! - `covalent-radii: pyykko-covalent`, the radii bond inference and the
//!   bond-length checks use;
//! - `vdw-radii: bondi`, the set the default analysis policy names, which
//!   surface area and clash scoring use;
//! - `metals: macromolecular-metals`, the sites surface area's default metal
//!   fallback reads.
//!
//! Standard-component annotation, selection and secondary-structure assignment
//! read no table. Everything else (other radius sets, ionic radii, validation
//! references, substitution matrices, component dictionaries) is not compiled
//! in: it is a choice a caller makes, and the profile or a registration makes it.
//!
//! The bytes are the repository's `data/` tables, copied beside this crate so a
//! published package carries them; a test keeps every copy identical to its
//! source. They decode through the same path as a cache file, once per process
//! on first selection, so every target that selects them gets the same tables.
//!
//! Without the feature the layer is absent and nothing is compiled in.

use crate::table_memory::MemoryTables;
use std::sync::Arc;

#[cfg(feature = "embedded-tables")]
const PROFILE: &str = "bond-reach: macromolecular-reach
covalent-radii: pyykko-covalent
vdw-radii: bondi
metals: macromolecular-metals
";

#[cfg(feature = "embedded-tables")]
const ENTRIES: [(&str, &str, &[u8]); 4] = [
    (
        "macromolecular-reach",
        <crate::BondReach as crate::TableKind>::EXTENSION,
        include_bytes!("../../data/macromolecular-reach.reach"),
    ),
    (
        "pyykko-covalent",
        <crate::CovalentRadii as crate::TableKind>::EXTENSION,
        include_bytes!("../../data/pyykko-covalent.covalent"),
    ),
    (
        "bondi",
        <crate::VdwRadii as crate::TableKind>::EXTENSION,
        include_bytes!("../../data/bondi.radii"),
    ),
    (
        "macromolecular-metals",
        <crate::MetalSites as crate::TableKind>::EXTENSION,
        include_bytes!("../../data/macromolecular-metals.metals"),
    ),
];

/// The compiled-in layer, built once per process.
#[cfg(feature = "embedded-tables")]
pub(crate) fn bundled() -> Option<Arc<MemoryTables>> {
    static TABLES: std::sync::OnceLock<Option<Arc<MemoryTables>>> = std::sync::OnceLock::new();
    TABLES
        .get_or_init(|| {
            let mut tables = MemoryTables::new("embedded");
            tables.set_profile(PROFILE);
            for (name, extension, bytes) in ENTRIES {
                tables.insert(name, extension, crate::TableBytes::from_static(bytes));
            }
            Some(Arc::new(tables))
        })
        .clone()
}

/// No layer: the library compiles in no table without the feature.
#[cfg(not(feature = "embedded-tables"))]
pub(crate) const fn bundled() -> Option<Arc<MemoryTables>> {
    None
}

#[cfg(test)]
#[path = "table_bundled_tests.rs"]
mod tests;

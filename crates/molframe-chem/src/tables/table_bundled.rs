//! The tables compiled into the library, when the `embedded-tables` feature asks.
//!
//! A read perceives bonds only when a bond-reach table is selected, and a
//! target without a filesystem (a browser) or a clean machine has no profile
//! to select one from. With the feature, the store's last layer is this set: the
//! minimum the default read-time perception needs, which is the bonding reach
//! table alone, and a profile naming it. Standard-component annotation and
//! secondary-structure assignment read no table.
//!
//! The bytes are the repository's `data/macromolecular-reach.reach`, copied
//! beside this crate so a published package carries them; a test keeps the two
//! identical. They decode through the same path as the cache file, so every
//! target that selects them perceives the same bonds.
//!
//! Without the feature the layer is absent and nothing is compiled in.

use crate::table_memory::MemoryTables;
use std::sync::Arc;

#[cfg(feature = "embedded-tables")]
const PROFILE: &str = "bond-reach: macromolecular-reach\n";

#[cfg(feature = "embedded-tables")]
const ENTRIES: [(&str, &str, &[u8]); 1] = [(
    "macromolecular-reach",
    <crate::BondReach as crate::TableKind>::EXTENSION,
    include_bytes!("../../data/macromolecular-reach.reach"),
)];

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

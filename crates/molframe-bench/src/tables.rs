//! A table store over the repository's own `data/` directory, for tests and benches.
//!
//! The library ships no scientific table, so every computation that needs one
//! is given a store. This one reads the files in `data/` and a profile that
//! selects the radii, reach, metal and Shannon tables.
//!
//! Cost: the store is built once per process; each call after is a clone of
//! a static reference.

use molframe_chem::{BondReach, Provided, TableRequest, TableStore};
use std::path::PathBuf;
use std::sync::{Arc, OnceLock};

/// The store over `data/` with `data/profile.txt` of this crate as its profile.
#[must_use]
pub fn table_store() -> &'static TableStore {
    shared_table_store()
}

/// The same store, shared.
#[must_use]
pub fn shared_table_store() -> &'static Arc<TableStore> {
    static STORE: OnceLock<Arc<TableStore>> = OnceLock::new();
    STORE.get_or_init(|| {
        let data = PathBuf::from(format!("{}/../../data", env!("CARGO_MANIFEST_DIR")));
        let profile = PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/data/profile.txt"));
        Arc::new(TableStore::new(data).with_profile(profile))
    })
}

/// The bonding reach table the profile selects, decoded once.
///
/// # Panics
///
/// When the profile's table is refused, which no test can recover from.
#[must_use]
pub fn bond_reach() -> &'static BondReach {
    static REACH: OnceLock<Arc<Provided<BondReach>>> = OnceLock::new();
    let provided = REACH.get_or_init(|| {
        match TableRequest::<BondReach>::profile().resolve(table_store(), "--reach") {
            Ok(provided) => provided,
            Err(error) => panic!("bond reach: {error}"),
        }
    });
    &provided.table
}

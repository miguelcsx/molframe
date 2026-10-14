//! Test helpers: decode the tables shipped in the repository's `data/` directory.

use crate::{
    BondReach, IonicRadii, MetalSites, MetalTables, Provided, RadiusSource, TableIdentity,
    TableKind, VdwRadii,
};
use std::sync::{Arc, OnceLock};

/// Decodes `data/<name>`, panicking with the reason when it is refused.
pub(crate) fn shipped<T: TableKind>(name: &str) -> Provided<T> {
    let path = format!("{}/../../data/{name}", env!("CARGO_MANIFEST_DIR"));
    let bytes = match std::fs::read(&path) {
        Ok(bytes) => bytes,
        Err(error) => panic!("{path}: {error}"),
    };
    match T::decode(&bytes, TableIdentity::supplied("x", "x", "x", "x")) {
        Ok(provided) => provided,
        Err(error) => panic!("{path}: {error}"),
    }
}

/// The shipped Bondi table as a radius source.
pub(crate) fn bondi() -> RadiusSource {
    RadiusSource::from(shipped::<VdwRadii>("bondi.radii"))
}

/// The shipped bonding reach table, decoded once.
pub(crate) fn reach() -> &'static BondReach {
    static REACH: OnceLock<BondReach> = OnceLock::new();
    REACH.get_or_init(|| shipped::<BondReach>("macromolecular-reach.reach").table)
}

/// The shipped metal and Shannon tables.
pub(crate) fn metal_tables() -> MetalTables {
    MetalTables {
        sites: Arc::new(shipped::<MetalSites>("macromolecular-metals.metals")),
        ionic: Some(Arc::new(shipped::<IonicRadii>("shannon-ionic.ionic"))),
    }
}

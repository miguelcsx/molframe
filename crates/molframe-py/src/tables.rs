//! Reference-table resolution for the bindings: one store, one refusal.
//!
//! Every binding that needs a radius or reach table resolves it here, from the
//! caller's `radii` argument (a cache entry name or a file path) when given and
//! from the profile otherwise. A table that cannot be supplied raises the
//! typed `E7102` or `E2201` exception, which names how to supply it.
//!
//! Cost: one profile stat and, on first use, one table decode per call.

use molframe::chemistry::{
    BondReach, CovalentRadii, Provided, RadiusSource, TableRequest, TableStore, VdwRadii,
};
use pyo3::PyResult;
use std::sync::Arc;

/// The van der Waals radii a call measures with.
pub(crate) fn radius_source(radii: Option<&str>) -> PyResult<RadiusSource> {
    RadiusSource::resolve(
        &TableRequest::<VdwRadii>::choose(radii),
        &TableStore::global(),
    )
    .map_err(crate::error::kernel)
}

/// The covalent radii the profile selects.
pub(crate) fn covalent_radii() -> PyResult<Arc<Provided<CovalentRadii>>> {
    TableRequest::<CovalentRadii>::profile()
        .resolve(&TableStore::global(), "covalent radii")
        .map_err(crate::error::kernel)
}

/// The bond-reach table the profile selects.
pub(crate) fn bond_reach() -> PyResult<Arc<Provided<BondReach>>> {
    TableRequest::<BondReach>::profile()
        .resolve(&TableStore::global(), "bond reach")
        .map_err(crate::error::kernel)
}

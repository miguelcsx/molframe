//! Reference-table resolution for the bindings: one store, one refusal.
//!
//! Every binding that needs a radius or reach table resolves it here, from the
//! caller's `radii` argument (a cache entry name or a file path) when given, from the
//! policy's `vdw_radii` when a policy is, and from the profile otherwise. A table that cannot be supplied raises the
//! typed `E7102` or `E2201` exception, which names how to supply it.
//!
//! Cost: one profile stat and, on first use, one table decode per call.

use molframe::chemistry::{
    BondReach, CovalentRadii, Provided, RadiusSelection, RadiusSource, TableRequest, TableStore,
    VdwRadii,
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

/// The van der Waals radii a call measures with when it may also carry a policy.
///
/// The one rule of [`RadiusSelection`], applied in Rust: the `radii` argument, else the
/// table the policy's `vdw_radii` names, else (no policy) the profile; an argument that
/// names another table than the policy is the typed `E7103` error. Returns the selection
/// so the kernel can check the table it runs with.
pub(crate) fn radii_under(
    radii: Option<&str>,
    policy: Option<&crate::policy::PyAnalysisPolicy>,
) -> PyResult<(RadiusSelection, RadiusSource)> {
    let selection = RadiusSelection::choose(radii, policy.is_some());
    let set = match policy {
        Some(policy) => policy.0.vdw_radii,
        None => molframe::AnalysisPolicy::default().vdw_radii,
    };
    let source = selection
        .resolve(set, &TableStore::global())
        .map_err(crate::error::kernel)?;
    Ok((selection, source))
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

//! Single checked-in operation catalog used by registration and parity checks.

use pyo3::prelude::*;

mod data;
mod structural;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Capability {
    pub(crate) name: &'static str,
    pub(crate) domain: &'static str,
    pub(crate) feature: &'static str,
    pub(crate) inputs: &'static str,
    pub(crate) result: &'static str,
    pub(crate) policy: bool,
    pub(crate) eager: bool,
    pub(crate) workflow: bool,
    pub(crate) execution_needs: &'static str,
    pub(crate) cost: molframe::Cost,
}

/// Every catalogued operation, grouped by what it operates on.
const GROUPS: [&[Capability]; 2] = [structural::STRUCTURAL, data::DATA];

/// Every catalogued operation.
pub(crate) fn capabilities() -> impl Iterator<Item = &'static Capability> {
    GROUPS.into_iter().flatten()
}

pub(crate) fn validate_registration(module: &Bound<'_, PyModule>) -> PyResult<()> {
    for capability in capabilities() {
        debug_assert!(!capability.inputs.is_empty());
        debug_assert!(!capability.result.is_empty());
        debug_assert!(!capability.feature.is_empty());
        debug_assert!(!capability.execution_needs.is_empty());
        debug_assert!(capability.eager || capability.workflow);
        let domain = module.getattr(capability.domain)?;
        if capability.eager {
            domain.getattr(capability.name)?;
        }
        let _ = (capability.policy, capability.cost);
    }
    Ok(())
}

//! Single checked-in operation catalog used by registration and parity checks.

use pyo3::prelude::*;

mod analysis;
mod audit;
mod chemistry;
mod compare;
mod crystal;
mod data;
mod formats;
mod geometry;
mod interop;
mod motif;
mod sequence;
mod spatial;
mod surface;
mod trajectory;
mod validation;

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

/// The namespace a block of catalogue rows belongs to.
#[derive(Clone, Copy)]
struct Domain(&'static str);

impl Domain {
    /// An eager, materialising, non-policy, non-workflow operation gated by the feature named
    /// after its domain.
    const fn row(
        self,
        name: &'static str,
        inputs: &'static str,
        result: &'static str,
        execution_needs: &'static str,
    ) -> Capability {
        Capability {
            name,
            domain: self.0,
            feature: self.0,
            inputs,
            result,
            policy: false,
            eager: true,
            workflow: false,
            execution_needs,
            cost: molframe::Cost::Materialize,
        }
    }
}

impl Capability {
    const fn cost(mut self, cost: molframe::Cost) -> Self {
        self.cost = cost;
        self
    }

    const fn feature(mut self, feature: &'static str) -> Self {
        self.feature = feature;
        self
    }

    const fn policy(mut self) -> Self {
        self.policy = true;
        self
    }

    const fn workflow(mut self) -> Self {
        self.workflow = true;
        self
    }
}

/// Every catalogued operation, grouped by what it operates on.
const GROUPS: [&[Capability]; 15] = [
    analysis::ROWS,
    audit::ROWS,
    chemistry::ROWS,
    compare::ROWS,
    crystal::ROWS,
    data::ROWS,
    formats::ROWS,
    geometry::ROWS,
    interop::ROWS,
    motif::ROWS,
    sequence::ROWS,
    spatial::ROWS,
    surface::ROWS,
    trajectory::ROWS,
    validation::ROWS,
];

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

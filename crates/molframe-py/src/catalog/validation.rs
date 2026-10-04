//! Catalogued operations of the `validation` namespace.

use super::Capability;

pub(super) const ROWS: &[Capability] = &[Capability {
    name: "clashes",
    domain: "validation",
    feature: "validation",
    inputs: "structure",
    result: "ClashTable",
    policy: false,
    eager: true,
    workflow: false,
    execution_needs: "cpu",
    cost: molframe::Cost::Materialize,
}];

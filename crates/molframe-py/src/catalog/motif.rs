//! Catalogued operations of the `motif` namespace.

use super::Capability;

pub(super) const ROWS: &[Capability] = &[Capability {
    name: "evaluate",
    domain: "motif",
    feature: "motif",
    inputs: "structure,specification,components,version,limits",
    result: "MotifReport",
    policy: true,
    eager: true,
    workflow: false,
    execution_needs: "cpu,memory",
    cost: molframe::Cost::Materialize,
}];

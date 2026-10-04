//! Catalogued operations of the `spatial` namespace.

use super::Capability;

pub(super) const ROWS: &[Capability] = &[Capability {
    name: "neighbor_pairs",
    domain: "spatial",
    feature: "spatial",
    inputs: "coordinates",
    result: "arrays",
    policy: false,
    eager: true,
    workflow: false,
    execution_needs: "cpu,spatial,memory",
    cost: molframe::Cost::Materialize,
}];

//! Catalogued operations of the `surface` namespace.

use super::Capability;

pub(super) const ROWS: &[Capability] = &[
    Capability {
        name: "sasa",
        domain: "surface",
        feature: "surface",
        inputs: "coordinates, radii",
        result: "array",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "lee_richards",
        domain: "surface",
        feature: "surface",
        inputs: "coordinates, radii",
        result: "array",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "cavities",
        domain: "surface",
        feature: "surface",
        inputs: "coordinates, radii",
        result: "list",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
];

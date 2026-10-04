//! Catalogued operations of the `compare` namespace.

use super::Capability;

pub(super) const ROWS: &[Capability] = &[
    Capability {
        name: "tm_score",
        domain: "compare",
        feature: "compare",
        inputs: "coordinates, coordinates",
        result: "float",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "gdt_ts",
        domain: "compare",
        feature: "compare",
        inputs: "coordinates, coordinates",
        result: "float",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "gdt_ha",
        domain: "compare",
        feature: "compare",
        inputs: "coordinates, coordinates",
        result: "float",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "weighted_rmsd",
        domain: "compare",
        feature: "compare",
        inputs: "coordinates, coordinates, weights",
        result: "float",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
];

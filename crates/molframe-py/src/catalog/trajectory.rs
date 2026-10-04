//! Catalogued operations of the `trajectory` namespace.

use super::Capability;

pub(super) const ROWS: &[Capability] = &[
    Capability {
        name: "read",
        domain: "trajectory",
        feature: "trajectory",
        inputs: "path",
        result: "Trajectory",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu,memory",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "rmsd",
        domain: "trajectory",
        feature: "trajectory",
        inputs: "positions",
        result: "Analysis",
        policy: true,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
];

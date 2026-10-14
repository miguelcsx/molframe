//! Catalogued operations of the `data` namespace.

use super::Capability;

pub(super) const ROWS: &[Capability] = &[
    Capability {
        name: "list",
        domain: "data",
        feature: "data",
        inputs: "manifest,directory?",
        result: "list[DataEntry]",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "verify",
        domain: "data",
        feature: "data",
        inputs: "manifest,directory?,name?",
        result: "list[DataEntry]",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "fetch",
        domain: "data",
        feature: "data",
        inputs: "manifest,name,limits,directory?,replace?",
        result: "Fetched",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu,network",
        cost: molframe::Cost::Materialize,
    },
];

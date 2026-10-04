//! Catalogued operations of the `formats` namespace.

use super::Capability;

pub(super) const ROWS: &[Capability] = &[
    Capability {
        name: "to_mmcif",
        domain: "formats",
        feature: "formats",
        inputs: "structure",
        result: "text",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu,memory",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "to_pdb",
        domain: "formats",
        feature: "formats",
        inputs: "structure",
        result: "text",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu,memory",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "to_bcif",
        domain: "formats",
        feature: "formats",
        inputs: "structure",
        result: "bytes",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu,memory",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "write",
        domain: "formats",
        feature: "formats",
        inputs: "structure, path",
        result: "file",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu,memory",
        cost: molframe::Cost::Materialize,
    },
];

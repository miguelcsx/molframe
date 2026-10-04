//! Catalogued operations of the `interop` namespace.

use super::Capability;

pub(super) const ROWS: &[Capability] = &[
    Capability {
        name: "coordinates",
        domain: "interop",
        feature: "interop",
        inputs: "structure",
        result: "CoordinateTensor",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "graph",
        domain: "interop",
        feature: "interop",
        inputs: "structure,nodes,edges,direction,cutoff,neighbors,node_features,edge_features,missing,backend,periodic",
        result: "Graph",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu,spatial,memory",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "write_atoms_ipc",
        domain: "interop",
        feature: "interop",
        inputs: "structure,path,metadata",
        result: "None",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu,io",
        cost: molframe::Cost::Materialize,
    },
    Capability {
        name: "write_atoms_parquet",
        domain: "interop",
        feature: "interop",
        inputs: "structure,path,metadata",
        result: "None",
        policy: false,
        eager: true,
        workflow: false,
        execution_needs: "cpu,io",
        cost: molframe::Cost::Materialize,
    },
];

//! Catalogued operations of the `interop` namespace.

use super::{Capability, Domain};

const D: Domain = Domain("interop");

pub(super) const ROWS: &[Capability] = &[
    D.row("coordinates", "structure", "CoordinateTensor", "cpu"),
    D.row("graph", "structure,nodes,edges,direction,cutoff,neighbors,node_features,edge_features,missing,backend,periodic", "Graph", "cpu,spatial,memory"),
    D.row("write_atoms_ipc", "structure,path,metadata", "None", "cpu,io"),
    D.row("write_atoms_parquet", "structure,path,metadata", "None", "cpu,io"),
];

//! Catalogued operations of the `spatial` namespace.

use super::{Capability, Domain};

const D: Domain = Domain("spatial");

pub(super) const ROWS: &[Capability] = &[
    D.row(
        "cross_pairs",
        "coordinates,coordinates",
        "arrays",
        "cpu,spatial,memory",
    ),
    D.row("knn", "coordinates", "arrays", "cpu,spatial,memory"),
    D.row(
        "neighbor_pairs",
        "coordinates",
        "arrays",
        "cpu,spatial,memory",
    ),
    D.row(
        "neighbor_pairs_chunks",
        "coordinates",
        "arrays",
        "cpu,spatial,memory",
    ),
];

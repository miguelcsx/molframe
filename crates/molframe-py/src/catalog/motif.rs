//! Catalogued operations of the `motif` namespace.

use super::{Capability, Domain};
use molframe::Cost;

const D: Domain = Domain("motif");

pub(super) const ROWS: &[Capability] = &[
    D.row(
        "evaluate",
        "structure,specification,components,version,limits",
        "MotifReport",
        "cpu,memory",
    )
    .policy(),
    D.row(
        "count",
        "structure,specification,components,version",
        "int",
        "cpu,memory",
    )
    .cost(Cost::Borrow)
    .policy(),
    D.row(
        "matches",
        "structure,specification,components,version",
        "MotifMatches",
        "cpu,memory",
    )
    .policy(),
    D.row(
        "first",
        "structure,specification,components,version",
        "MotifMatches",
        "cpu,memory",
    )
    .policy(),
    D.row(
        "iter_matches",
        "structure,specification,components,version",
        "MatchIterator",
        "cpu,memory",
    )
    .policy(),
];

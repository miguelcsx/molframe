//! Catalogued operations of the `data` namespace.

use super::{Capability, Domain};

const D: Domain = Domain("data");

pub(super) const ROWS: &[Capability] = &[
    D.row("list", "manifest,directory?", "list[DataEntry]", "cpu"),
    D.row(
        "verify",
        "manifest,directory?,name?",
        "list[DataEntry]",
        "cpu",
    ),
    D.row(
        "fetch",
        "manifest,name,limits,directory?,replace?",
        "Fetched",
        "cpu,network",
    ),
];

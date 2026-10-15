//! Catalogued operations of the `audit` namespace.

use super::{Capability, Domain};

const D: Domain = Domain("audit");

pub(super) const ROWS: &[Capability] = &[D
    .row(
        "run",
        "plan,analyse,metric,project",
        "AuditResult",
        "cpu,memory",
    )
    .policy()];

//! Catalogued operations of the `audit` namespace.

use super::Capability;

pub(super) const ROWS: &[Capability] = &[Capability {
    name: "run",
    domain: "audit",
    feature: "audit",
    inputs: "plan,analyse,metric,project",
    result: "AuditResult",
    policy: true,
    eager: true,
    workflow: false,
    execution_needs: "cpu,memory",
    cost: molframe::Cost::Materialize,
}];

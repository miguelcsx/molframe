//! The complete table of declared cases.

use crate::harness::Case;

pub fn all() -> Vec<Case> {
    let mut cases = Vec::new();
    cases.extend(crate::spatial_cases::cases());
    cases.extend(crate::geom_cases::cases());
    cases.extend(crate::compare_cases::cases());
    cases.extend(crate::traj_cases::cases());
    cases
}

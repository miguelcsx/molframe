//! Stable diagnostic codes for the ways an audit plan can be refused.

use crate::PlanError;
use molframe_core::{Code, diagnostic_from};

diagnostic_from!(PlanError, |error| match error {
    PlanError::EmptyDimension(_) | PlanError::DuplicateDimension(_) => Code::E5101,
    PlanError::Conflict { .. } => Code::E6004,
    PlanError::NotRead(_) => Code::E6103,
    PlanError::CostOverflow | PlanError::LimitExceeded { .. } => Code::E7001,
});

//! Shared support for structure operation execution.

use super::plan::value::ExecutionPlanError;
use molframe_core::structure::Structure;
use molframe_spatial::PeriodicBox;

pub(super) fn periodic_box(
    structure: &Structure,
    periodic: bool,
) -> Result<Option<PeriodicBox>, ExecutionPlanError> {
    if !periodic {
        return Ok(None);
    }
    let cell = structure
        .data()
        .cell
        .ok_or_else(|| ExecutionPlanError::Governed("periodic GNM requires a unit cell".into()))?;
    PeriodicBox::from_cell(cell)
        .map(Some)
        .map_err(|error| ExecutionPlanError::Governed(error.to_string().into()))
}

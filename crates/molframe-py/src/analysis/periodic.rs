//! Shared conversion from a structure cell to explicit periodic geometry.

pub(super) fn periodic_box(
    structure: &molframe::Structure,
    periodic: bool,
) -> Result<Option<molframe::PeriodicBox>, String> {
    if !periodic {
        return Ok(None);
    }
    let cell = structure
        .data()
        .cell
        .ok_or_else(|| "periodic analysis requires a unit cell".to_owned())?;
    molframe::PeriodicBox::from_cell(cell)
        .map(Some)
        .map_err(|error| error.to_string())
}

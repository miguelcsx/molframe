//! Stable trajectory summaries shared by human and machine output.

use pdbiox::traj::{TrajectoryData, TrajectoryFormat, Units};

#[derive(Debug, PartialEq)]
pub(super) struct Summary {
    pub format: &'static str,
    pub frames: usize,
    pub atoms: usize,
    pub first_time: Option<f64>,
    pub last_time: Option<f64>,
    pub frames_with_cell: usize,
    pub frames_with_velocities: usize,
    pub frames_with_forces: usize,
}

impl Summary {
    pub(super) fn from_data(data: &TrajectoryData) -> Self {
        Self {
            format: format_name(data.format),
            frames: data.frames.len(),
            atoms: data.frames.first().map_or(0, |frame| frame.positions.len()),
            first_time: data.frames.first().and_then(|frame| frame.time),
            last_time: data.frames.last().and_then(|frame| frame.time),
            frames_with_cell: count(data, |frame| frame.cell.is_some()),
            frames_with_velocities: count(data, |frame| frame.velocities.is_some()),
            frames_with_forces: count(data, |frame| frame.forces.is_some()),
        }
    }
}

fn count(data: &TrajectoryData, predicate: impl Fn(&pdbiox::traj::Timestep) -> bool) -> usize {
    data.frames.iter().filter(|frame| predicate(frame)).count()
}

pub(super) const fn format_name(format: TrajectoryFormat) -> &'static str {
    match format {
        TrajectoryFormat::Xtc => "xtc",
        TrajectoryFormat::Trr => "trr",
        TrajectoryFormat::Dcd => "dcd",
        TrajectoryFormat::AmberNetcdf => "amber-netcdf",
        TrajectoryFormat::Tng => "tng",
        TrajectoryFormat::Gsd => "gsd",
        TrajectoryFormat::H5md => "h5md",
        TrajectoryFormat::Trz => "trz",
        _ => "unsupported",
    }
}

pub(super) const fn canonical_length_unit() -> &'static str {
    Units::CANONICAL.length
}

pub(super) const fn canonical_time_unit() -> &'static str {
    Units::CANONICAL.time
}

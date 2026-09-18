//! Process-isolated XTC pull-reader resource case.

use std::hint::black_box;
use std::path::Path;

use super::{ResourceRecord, measure_case};

pub(super) fn read_file(path: &Path) -> Result<ResourceRecord, String> {
    let mut model_count = 0_u64;
    let mut record = measure_case("xtc_file", || {
        let mut reader = molframe::traj::read_trajectory(
            path,
            &molframe::traj::TrajectoryReaderOptions::default(),
        )
        .map_err(|error| format!("{} XTC open failed: {error}", path.display()))?;
        let atoms = reader.n_atoms();
        let mut timestep = molframe::traj::Timestep::default();
        let mut frames = 0_u64;
        while reader
            .read_next(&mut timestep)
            .map_err(|error| format!("{} XTC read failed: {error}", path.display()))?
        {
            if timestep.positions.len() != atoms {
                return Err("XTC pull reader changed atom count".to_owned());
            }
            black_box(&timestep.positions);
            frames = frames
                .checked_add(1)
                .ok_or_else(|| "XTC frame count exceeds u64".to_owned())?;
        }
        model_count = frames;
        let atoms = u64::try_from(atoms).map_err(|_| "XTC atom count exceeds u64".to_owned())?;
        frames
            .checked_mul(1_000_000_000)
            .and_then(|digest| digest.checked_add(atoms))
            .ok_or_else(|| "XTC benchmark digest exceeds u64".to_owned())
    })?;
    record.model_count = Some(model_count);
    Ok(record)
}

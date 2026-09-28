//! Process-isolated whole-trajectory TNG resource cases.

use std::hint::black_box;
use std::path::Path;

use molframe::trajectory::{Timestep, TngWriteOptions, write_tng};
use molframe_bench::{Sample, coordinates, structure};
use tempfile::NamedTempFile;

use super::{ResourceRecord, measure_retained_case};

pub(super) fn synthetic() -> Result<ResourceRecord, String> {
    let positions = coordinates(&structure(Sample::Small));
    let frames: Vec<_> = (0..128)
        .map(|frame| Timestep {
            frame,
            time: Some(f64::from(
                u32::try_from(frame).expect("synthetic frame index fits in u32"),
            )),
            positions: positions.clone(),
            ..Timestep::default()
        })
        .collect();
    let file =
        NamedTempFile::new().map_err(|error| format!("TNG temporary file failed: {error}"))?;
    write_tng(file.path(), &frames, TngWriteOptions::default())
        .map_err(|error| format!("TNG fixture write failed: {error}"))?;
    read("tng_synthetic_128x660", file.path())
}

pub(super) fn read_file(path: &Path) -> Result<ResourceRecord, String> {
    read("tng_file", path)
}

fn read(name: &'static str, path: &Path) -> Result<ResourceRecord, String> {
    measure_retained_case(name, || {
        let trajectory = molframe::trajectory::parse_tng(path)
            .map_err(|error| format!("{} TNG read failed: {error}", path.display()))?;
        let frames = u64::try_from(trajectory.frames.len())
            .map_err(|_| "TNG frame count exceeds u64".to_owned())?;
        let atoms = match trajectory.frames.first() {
            Some(frame) => u64::try_from(frame.positions.len())
                .map_err(|_| "TNG atom count exceeds u64".to_owned())?,
            None => 0,
        };
        let digest = frames
            .checked_mul(1_000_000_000)
            .and_then(|digest| digest.checked_add(atoms))
            .ok_or_else(|| "TNG benchmark digest exceeds u64".to_owned())?;
        black_box(&trajectory);
        Ok((digest, Some(frames), trajectory))
    })
}

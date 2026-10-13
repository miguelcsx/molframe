//! Trajectory file reads: frame blocks decoded by one worker or several.
//!
//! The Huge tier is 100 000 frames. A case at one worker is the serial
//! baseline; the same case at more workers shows what the shared frame index
//! and the interleaved block schedule buy. Checksums fold every block in file
//! order, so rows at different worker counts must agree bit for bit.

use std::path::{Path, PathBuf};

use crate::harness::{Axis, Case, Outcome, Run};
use molframe_bench::{Checksum, Tier, ensemble};
use molframe_traj::{
    DcdReader, DcdWriteOptions, ParallelChunks, Timestep, TrajectoryError, TrrReader,
    TrrWriteOptions, XtcReader, XtcWriteOptions, write_dcd, write_trr, write_xtc,
};

const ATOMS: usize = 100;
const CHUNK_FRAMES: usize = 256;
const MEMORY_LIMIT: usize = 1 << 30;
const FILE_TIERS: &[Tier] = &[Tier::Small, Tier::Medium, Tier::Large, Tier::Huge];

pub fn cases() -> Vec<Case> {
    vec![
        Case {
            id: "traj.file_chunks.xtc.100_atoms",
            axis: Axis::Frames,
            tiers: FILE_TIERS,
            run: xtc_case,
        },
        Case {
            id: "traj.file_chunks.trr.100_atoms",
            axis: Axis::Frames,
            tiers: FILE_TIERS,
            run: trr_case,
        },
        Case {
            id: "traj.file_chunks.dcd.100_atoms",
            axis: Axis::Frames,
            tiers: FILE_TIERS,
            run: dcd_case,
        },
    ]
}

/// A scratch file removed when the case ends.
struct Scratch(PathBuf);

impl Drop for Scratch {
    fn drop(&mut self) {
        // The file is a cache of generated data; a failed removal costs only disk.
        drop(std::fs::remove_file(&self.0));
    }
}

#[allow(clippy::cast_precision_loss)]
fn timesteps(frames: usize) -> Vec<Timestep> {
    let ens = ensemble(ATOMS, frames, 0.8, 5);
    (0..frames)
        .map(|frame| Timestep {
            frame,
            time: Some(frame as f64 * 0.5),
            positions: ens.frame(frame).to_vec(),
            ..Timestep::default()
        })
        .collect()
}

fn scratch(name: &str, frames: usize, bytes: Vec<u8>) -> Result<Scratch, String> {
    let path = std::env::temp_dir().join(format!(
        "molframe-scale-{}-{name}-{frames}",
        std::process::id()
    ));
    std::fs::write(&path, bytes).map_err(|e| e.to_string())?;
    Ok(Scratch(path))
}

fn xtc_case(run: &Run<'_>, frames: usize) -> Result<Outcome, String> {
    let bytes =
        write_xtc(&timesteps(frames), XtcWriteOptions::default()).map_err(|e| e.to_string())?;
    let file = scratch("xtc.xtc", frames, bytes)?;
    XtcReader::open(&file.0).map_err(|e| e.to_string())?;
    chunk_case(run, frames, &file.0, |path| {
        XtcReader::open(path).map_err(|_| TrajectoryError::RandomAccessUnavailable)
    })
}

fn trr_case(run: &Run<'_>, frames: usize) -> Result<Outcome, String> {
    let bytes =
        write_trr(&timesteps(frames), TrrWriteOptions::default()).map_err(|e| e.to_string())?;
    let file = scratch("trr.trr", frames, bytes)?;
    TrrReader::open(&file.0, MEMORY_LIMIT).map_err(|e| e.to_string())?;
    chunk_case(run, frames, &file.0, |path| {
        TrrReader::open(path, MEMORY_LIMIT)
    })
}

fn dcd_case(run: &Run<'_>, frames: usize) -> Result<Outcome, String> {
    let bytes =
        write_dcd(&timesteps(frames), &DcdWriteOptions::default()).map_err(|e| e.to_string())?;
    let file = scratch("dcd.dcd", frames, bytes)?;
    DcdReader::open(&file.0).map_err(|e| e.to_string())?;
    chunk_case(run, frames, &file.0, |path| {
        DcdReader::open(path).map_err(|_| TrajectoryError::RandomAccessUnavailable)
    })
}

/// Reads every block of `path` through the parallel chunk iterator.
fn chunk_case<R>(
    run: &Run<'_>,
    frames: usize,
    path: &Path,
    open: impl Fn(&Path) -> Result<R, TrajectoryError>,
) -> Result<Outcome, String>
where
    R: molframe_traj::TrajectoryReader + Send + 'static,
{
    let sample = run.measure(
        || {
            let chunks = ParallelChunks::open(|| open(path), 0..frames, CHUNK_FRAMES, &run.context)
                .map_err(|e| e.to_string())?;
            let mut sum = Checksum::default();
            let mut seen = 0_usize;
            for chunk in chunks {
                let chunk = chunk.map_err(|e| e.to_string())?;
                seen += chunk.indices.len();
                let folded = chunk
                    .coordinates
                    .iter()
                    .fold(0_u64, |acc, v| acc.wrapping_add(u64::from(v.to_bits())));
                sum.word(folded);
                sum.word(chunk.indices.first().map_or(0, |first| *first));
            }
            if seen != frames {
                return Err(format!("read {seen} of {frames} frames"));
            }
            Ok(sum.value())
        },
        |value| *value,
    )?;
    Ok(Outcome {
        sample,
        items: (frames * ATOMS) as u64,
    })
}

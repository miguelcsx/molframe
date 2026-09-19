//! Decode throughput at a scale the 4x4 rows cannot show: 128 frames of the
//! Small fixture. The `Throughput::Bytes` rows report megabytes per second,
//! and the XTC and TNG rows are the evidence base for the buffer-reuse fixes
//! in those readers.

use crate::BenchRequired as _;
use crate::synthetic_frames;
use crate::usize_to_u32;
use criterion::{Criterion, Throughput, black_box};
use molframe_bench::{Sample, coordinates, structure};
use molframe_traj::{
    DcdWriteOptions, Timestep, TngWriteOptions, TrrWriteOptions, XtcWriteOptions, parse_dcd,
    parse_tng, parse_trr, parse_xtc, write_dcd, write_tng, write_trr, write_xtc,
};

/// Decode throughput at a scale the 4x4 rows cannot show: 128 frames of the
/// Small fixture. The `Throughput::Bytes` rows report megabytes per second,
/// and the XTC and TNG rows are the evidence base for the buffer-reuse fixes
/// in those readers.
pub(crate) fn bench_decode_throughput(c: &mut Criterion) {
    let sample = structure(Sample::Small);
    let frames: Vec<_> = synthetic_frames(&coordinates(&sample), 128)
        .into_iter()
        .enumerate()
        .map(|(frame, positions)| Timestep {
            frame,
            time: Some(f64::from(usize_to_u32(frame))),
            positions,
            ..Timestep::default()
        })
        .collect();

    let xtc = write_xtc(&frames, XtcWriteOptions::default()).required("XTC fixture failed");
    let trr = write_trr(&frames, TrrWriteOptions::default()).required("TRR fixture failed");
    let dcd = write_dcd(&frames, &DcdWriteOptions::default()).required("DCD fixture failed");
    let tng_file = tempfile::NamedTempFile::new().required("TNG temporary file failed");
    write_tng(tng_file.path(), &frames, TngWriteOptions::default()).required("TNG fixture failed");
    let tng_path = tng_file.path().to_path_buf();

    let mut group = c.benchmark_group("traj_decode");
    for (label, bytes) in [("xtc", &xtc), ("trr", &trr), ("dcd", &dcd)] {
        group.throughput(Throughput::Bytes(bytes.len() as u64));
        group.bench_function(format!("parse/{label}/128x{}", sample.atom_count()), |b| {
            b.iter(|| {
                let frame_count = match label {
                    "xtc" => parse_xtc(&xtc)
                        .map(|value| value.frames.len())
                        .map_err(|error| error.to_string()),
                    "trr" => parse_trr(&trr)
                        .map(|value| value.frames.len())
                        .map_err(|error| error.to_string()),
                    _ => parse_dcd(&dcd)
                        .map(|value| value.frames.len())
                        .map_err(|error| error.to_string()),
                };
                let count = match frame_count {
                    Ok(count) => count,
                    Err(error) => panic!("decode benchmark failed: {error}"),
                };
                black_box(count);
            });
        });
    }
    let tng_bytes = std::fs::read(&tng_path).required("TNG benchmark read failed");
    group.throughput(Throughput::Bytes(tng_bytes.len() as u64));
    group.bench_function(format!("parse/tng/128x{}", sample.atom_count()), |b| {
        b.iter(|| {
            let parsed = parse_tng(&tng_path).required("TNG decode benchmark failed");
            black_box(parsed.frames.len());
        });
    });
    group.finish();
}

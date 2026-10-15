//! Cartesian PCA at scales where the covariance route and the snapshot route
//! each dominate: frames much greater than 3N, and frames much smaller.

use criterion::{Criterion, black_box};
use molframe_traj::{CartesianFit, cartesian_pca};

const MEMORY: usize = 8 * 1024 * 1024 * 1024;

/// Deterministic low-rank-plus-noise frames, `frames x atoms`.
pub(crate) fn low_rank_frames(frames: usize, atoms: usize) -> Vec<Vec<[f32; 3]>> {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut next = move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        f32::from(top_sixteen_bits(state)) / 65_536.0 - 0.5
    };
    let base: Vec<[f32; 3]> = (0..atoms).map(|_| [next(), next(), next()]).collect();
    let modes: Vec<Vec<[f32; 3]>> = (0..4)
        .map(|_| (0..atoms).map(|_| [next(), next(), next()]).collect())
        .collect();
    (0..frames)
        .map(|_| {
            let weights: Vec<f32> = [4.0_f32, 3.0, 2.0, 1.0]
                .into_iter()
                .map(|scale| next() * scale)
                .collect();
            (0..atoms)
                .map(|atom| {
                    let mut point = base[atom];
                    for (weight, mode) in weights.iter().zip(&modes) {
                        for axis in 0..3 {
                            point[axis] += weight * mode[atom][axis];
                        }
                    }
                    point.map(|value| value + 0.01 * next())
                })
                .collect()
        })
        .collect()
}

pub(crate) fn bench_pca_scale(c: &mut Criterion) {
    let mut group = c.benchmark_group("traj_pca_scale");
    group.sample_size(10);
    for (frames, atoms) in [(4_000, 100), (2_000, 500), (300, 5_000)] {
        let data = low_rank_frames(frames, atoms);
        group.bench_function(format!("cartesian_pca/{frames}x{atoms}"), |b| {
            b.iter(|| black_box(cartesian_pca(&data, CartesianFit::None, 5, MEMORY)));
        });
    }
    group.finish();
}

/// The high sixteen bits of `state`, read as bytes so no conversion can fail.
fn top_sixteen_bits(state: u64) -> u16 {
    let bytes = state.to_le_bytes();
    u16::from_le_bytes([bytes[6], bytes[7]])
}

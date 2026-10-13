//! Deterministic synthetic coordinates at any size.

use super::rng::SplitMix64;

/// Number density of heavy atoms in folded protein, per cubic ångström.
const PROTEIN_DENSITY: f32 = 0.08;

/// A coordinate set plus the box that contains it.
#[derive(Clone, Debug)]
pub struct Cloud {
    /// Atom positions.
    pub positions: Vec<[f32; 3]>,
    /// Edge length of the cube the points were drawn from.
    pub edge: f32,
}

/// `count` points uniform in a cube sized for protein-like density.
#[must_use]
pub fn uniform_cloud(count: usize, seed: u64) -> Cloud {
    let mut rng = SplitMix64::new(seed);
    let edge =
        (crate::numeric::f64_to_f32(crate::numeric::usize_to_f64(count)) / PROTEIN_DENSITY).cbrt();
    let positions = (0..count)
        .map(|_| [rng.unit() * edge, rng.unit() * edge, rng.unit() * edge])
        .collect();
    Cloud { positions, edge }
}

/// A persistent random walk with 3.8 Å steps: a C-alpha trace, or, at higher
/// counts, a long polymer. Successive steps turn by a bounded angle, so the
/// chain neither collapses nor runs straight.
#[must_use]
pub fn chain(count: usize, seed: u64) -> Vec<[f32; 3]> {
    let mut rng = SplitMix64::new(seed);
    let mut heading = rng.direction();
    let mut at = [0.0_f32; 3];
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        out.push(at);
        let kick = rng.direction();
        let mut next = [
            heading[0] + 0.9 * kick[0],
            heading[1] + 0.9 * kick[1],
            heading[2] + 0.9 * kick[2],
        ];
        let norm = (next[0] * next[0] + next[1] * next[1] + next[2] * next[2]).sqrt();
        for axis in &mut next {
            *axis /= norm;
        }
        heading = next;
        for axis in 0..3 {
            at[axis] += 3.8 * heading[axis];
        }
    }
    out
}

/// `positions` rotated by a seeded rotation, translated by a seeded vector.
#[must_use]
pub fn rotated(positions: &[[f32; 3]], seed: u64) -> Vec<[f32; 3]> {
    let mut rng = SplitMix64::new(seed);
    let first = rng.direction();
    let mut second = rng.direction();
    let along = first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
    for axis in 0..3 {
        second[axis] -= along * first[axis];
    }
    let norm = (second[0] * second[0] + second[1] * second[1] + second[2] * second[2])
        .sqrt()
        .max(1.0e-6);
    for axis in &mut second {
        *axis /= norm;
    }
    let third = [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ];
    let shift = [rng.unit() * 50.0, rng.unit() * 50.0, rng.unit() * 50.0];
    positions
        .iter()
        .map(|p| {
            [
                first[0] * p[0] + second[0] * p[1] + third[0] * p[2] + shift[0],
                first[1] * p[0] + second[1] * p[1] + third[1] * p[2] + shift[1],
                first[2] * p[0] + second[2] * p[1] + third[2] * p[2] + shift[2],
            ]
        })
        .collect()
}

/// A trajectory in one contiguous frame-major buffer.
#[derive(Clone, Debug)]
pub struct Ensemble {
    /// Atoms per frame.
    pub atoms: usize,
    /// Number of frames.
    pub frames: usize,
    /// `frames * atoms` positions, frame-major.
    pub positions: Vec<[f32; 3]>,
}

impl Ensemble {
    /// One frame's positions.
    #[must_use]
    pub fn frame(&self, index: usize) -> &[[f32; 3]] {
        &self.positions[index * self.atoms..(index + 1) * self.atoms]
    }
}

/// A thermal ensemble: a chain with per-frame Gaussian noise of `sigma` ångström.
#[must_use]
pub fn ensemble(atoms: usize, frames: usize, sigma: f32, seed: u64) -> Ensemble {
    let base = chain(atoms, seed);
    let mut rng = SplitMix64::new(seed ^ 0xA5A5_5A5A);
    let mut positions = Vec::with_capacity(atoms * frames);
    for _ in 0..frames {
        for p in &base {
            positions.push([
                p[0] + sigma * rng.normal(),
                p[1] + sigma * rng.normal(),
                p[2] + sigma * rng.normal(),
            ]);
        }
    }
    Ensemble {
        atoms,
        frames,
        positions,
    }
}

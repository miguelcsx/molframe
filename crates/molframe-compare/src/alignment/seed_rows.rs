//! One running best per seed, held in two flat arrays.
//!
//! A search keeps the best score and transform each seed reached on every
//! tracked slot. Storing those as a row of one arena, rather than two vectors
//! per seed, removes the per-seed allocations that dominated small pairs.

use crate::seed_objective::BestRef;
use crate::seed_search::Schedule;
use molframe_geom::Rigid;
use rayon::prelude::*;

/// `slots` entries per seed, written once by the worker that runs the seed and
/// read after the join.
pub(crate) struct SeedRows {
    slots: usize,
    scores: Vec<f64>,
    transforms: Vec<Option<Rigid>>,
}

impl SeedRows {
    /// Runs `work` for every seed with a fresh row and a worker-local scratch
    /// from `make_scratch`, placed by `schedule`.
    ///
    /// Each row depends only on its seed, so the arena is identical for either
    /// schedule and any worker count.
    pub(crate) fn run<S, T: Copy + Sync>(
        seeds: &[T],
        slots: usize,
        schedule: Schedule,
        make_scratch: impl Fn() -> S + Sync + Send,
        work: impl Fn(&mut S, T, BestRef<'_>) + Sync + Send,
    ) -> Self {
        let mut rows = Self {
            slots,
            scores: vec![-1.0; seeds.len() * slots],
            transforms: vec![None; seeds.len() * slots],
        };
        let job = |scratch: &mut S, ((scores, transforms), &seed): ((&mut [f64], &mut [_]), _)| {
            work(scratch, seed, BestRef::reset(scores, transforms));
        };
        match schedule {
            Schedule::Pool => rows
                .scores
                .par_chunks_mut(slots)
                .zip(rows.transforms.par_chunks_mut(slots))
                .zip(seeds.par_iter())
                .for_each_init(&make_scratch, job),
            Schedule::Inline => {
                let mut scratch = make_scratch();
                rows.scores
                    .chunks_mut(slots)
                    .zip(rows.transforms.chunks_mut(slots))
                    .zip(seeds)
                    .for_each(|item| job(&mut scratch, item));
            }
        }
        rows
    }

    /// The score a seed result holds for a slot, matched by transform.
    pub(crate) fn score_of(&self, slot: usize, transform: Rigid) -> f64 {
        (0..self.count())
            .find(|&seed| self.transforms(seed)[slot] == Some(transform))
            .map_or(0.0, |seed| self.scores(seed)[slot])
    }

    pub(crate) fn count(&self) -> usize {
        self.scores.len() / self.slots
    }

    pub(crate) fn scores(&self, seed: usize) -> &[f64] {
        &self.scores[seed * self.slots..(seed + 1) * self.slots]
    }

    pub(crate) fn transforms(&self, seed: usize) -> &[Option<Rigid>] {
        &self.transforms[seed * self.slots..(seed + 1) * self.slots]
    }
}

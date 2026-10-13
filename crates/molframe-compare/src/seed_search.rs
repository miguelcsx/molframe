//! The seed search shared by TM-score and GDT.
//!
//! Both scores are a *maximum over rigid superpositions* of a sum over atoms, so
//! a single least-squares fit is only one candidate. This module searches for
//! better ones in the way LGA and TM-align do: it seeds a superposition from a
//! contiguous fragment of the correspondence list, then iterates "fit the set,
//! keep the atoms that are close, refit" until the set stops changing, scoring
//! every superposition it visits.
//!
//! # Seeds
//!
//! Fragment lengths start at the whole list and halve down to
//! [`MIN_SEED`] atoms; each length slides in steps of a quarter of its length,
//! widened so that no length has more than [`MAX_SEEDS_PER_LENGTH`] starts. The
//! seed count is therefore bounded by `MAX_SEEDS_PER_LENGTH * log2(n)`
//! whatever `n` is.
//!
//! # Bounded cost
//!
//! Seeds are refined on at most [`SAMPLE_CAP`] atoms (a stride subsample, the
//! whole list when it is shorter). Up to that length the search is exact. Above
//! it the seed search costs a fixed amount, and only the few best seed
//! transforms per tracked score are refined and scored on every atom, so the
//! total grows linearly in `n`. Every number returned is the score of an actual
//! transform on all atoms, so it never exceeds the true maximum.
//!
//! # Determinism
//!
//! Seeds run in parallel, but each seed's result is collected in seed order
//! and reduced sequentially with a strict "greater wins" rule, so the answer
//! does not depend on the number of workers or on scheduling.

use molframe_geom::{Rigid, superpose};
use rayon::prelude::*;

use crate::numeric::usize_to_f64;
pub(crate) use crate::seed_objective::Objective;
use crate::seed_objective::{Best, BestRef};
use crate::seed_rows::SeedRows;

/// The most atoms the seed search refines on.
pub(crate) const SAMPLE_CAP: usize = 1024;

/// Refinement rounds per seed and track before giving up on convergence.
const MAX_ITERATIONS: usize = 20;

/// The most start positions tried for any one seed length.
const MAX_SEEDS_PER_LENGTH: usize = 128;

/// The shortest seed fragment, in atoms.
const MIN_SEED: usize = 4;

/// How much a fitting set smaller than three atoms is widened per attempt.
const WIDEN_STEP: f64 = 0.5;

/// Refinement rounds on the full atom list for a sampled-search candidate.
const POLISH_ITERATIONS: usize = 10;

/// Seed transforms per tracked score that are carried to the full atom list.
const CANDIDATES_PER_SLOT: usize = 8;

/// Where the independent seeds of one search run.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Schedule {
    /// Across the shared rayon pool, one scratch per task.
    Pool,
    /// On the calling thread with one scratch, for callers that already
    /// parallelise over many searches.
    Inline,
}

/// Two aligned coordinate lists.
#[derive(Clone, Copy)]
struct Problem<'a> {
    model: &'a [[f32; 3]],
    reference: &'a [[f32; 3]],
}

/// Buffers reused across every round of one worker.
struct Scratch {
    distances: Vec<f64>,
    scores: Vec<f64>,
    current: Vec<usize>,
    next: Vec<usize>,
    moving: Vec<[f32; 3]>,
    fixed: Vec<[f32; 3]>,
}

impl Scratch {
    fn new(length: usize, slots: usize) -> Self {
        Self {
            distances: vec![0.0; length],
            scores: vec![0.0; slots],
            current: Vec::with_capacity(length),
            next: Vec::with_capacity(length),
            moving: Vec::with_capacity(length),
            fixed: Vec::with_capacity(length),
        }
    }

    /// Distances and tracked scores of every atom under `transform`.
    fn evaluate(&mut self, problem: Problem<'_>, objective: &Objective, transform: Rigid) {
        for ((slot, &point), &target) in self
            .distances
            .iter_mut()
            .zip(problem.model)
            .zip(problem.reference)
        {
            *slot = molframe_geom::distance(transform.apply(point), target);
        }
        objective.score(&self.distances, &mut self.scores);
    }

    /// Fits the atoms in `current`, or `None` when they cannot be fitted.
    fn fit_current(&mut self, problem: Problem<'_>) -> Option<Rigid> {
        self.moving.clear();
        self.fixed.clear();
        for &index in &self.current {
            self.moving.push(problem.model[index]);
            self.fixed.push(problem.reference[index]);
        }
        superpose(&self.moving, &self.fixed)
            .ok()
            .map(|fit| fit.transform)
    }

    /// Fills `next` with the atoms within `cutoff` of the last evaluated
    /// transform, widening until three qualify. Returns whether a set was found.
    fn select_next(&mut self, cutoff: f64) -> bool {
        let mut limit = cutoff;
        let mut ceiling = None;
        loop {
            self.next.clear();
            for (index, distance) in self.distances.iter().enumerate() {
                if *distance <= limit {
                    self.next.push(index);
                }
            }
            if self.next.len() >= 3 {
                return true;
            }
            let top = *ceiling
                .get_or_insert_with(|| self.distances.iter().copied().fold(0.0_f64, f64::max));
            if limit >= top {
                return false;
            }
            limit += WIDEN_STEP;
        }
    }

    /// Iterates fit, score, select from the set already in `current`.
    ///
    /// `first` stands in for the fit of the initial set when it is known.
    fn refine(
        &mut self,
        problem: Problem<'_>,
        objective: &Objective,
        cutoff: f64,
        mut first: Option<Rigid>,
        iterations: usize,
        best: &mut BestRef<'_>,
    ) {
        for _ in 0..iterations {
            let fitted = match first.take() {
                Some(transform) => Some(transform),
                None => self.fit_current(problem),
            };
            let Some(transform) = fitted else {
                break;
            };
            self.evaluate(problem, objective, transform);
            best.offer(&self.scores, transform);
            if !self.select_next(cutoff) || self.next == self.current {
                break;
            }
            std::mem::swap(&mut self.current, &mut self.next);
        }
    }
}

/// The `(start, length)` of every seed fragment for `length` atoms.
fn seeds(length: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut fragment = length;
    loop {
        let last_start = length - fragment;
        let step = (fragment / 4)
            .max(1)
            .max(last_start.div_ceil(MAX_SEEDS_PER_LENGTH));
        let mut start = 0;
        while start <= last_start {
            out.push((start, fragment));
            start += step;
        }
        // The last window is always tried, whatever the step did.
        if !last_start.is_multiple_of(step) {
            out.push((last_start, fragment));
        }
        if fragment / 2 < MIN_SEED {
            break;
        }
        fragment /= 2;
    }
    out
}

/// Refines one seed fragment on every track.
fn run_seed(
    problem: Problem<'_>,
    objective: &Objective,
    scratch: &mut Scratch,
    seed: (usize, usize),
    mut best: BestRef<'_>,
) {
    let (start, fragment) = seed;
    let range = start..start + fragment;
    let Ok(fit) = superpose(
        &problem.model[range.clone()],
        &problem.reference[range.clone()],
    ) else {
        return;
    };
    for &cutoff in &objective.select {
        scratch.current.clear();
        scratch.current.extend(range.clone());
        scratch.refine(
            problem,
            objective,
            cutoff,
            Some(fit.transform),
            MAX_ITERATIONS,
            &mut best,
        );
    }
}

/// The outcome of a search, scored on every atom.
#[derive(Debug)]
pub(crate) struct Outcome {
    /// The largest score found per tracked slot (counts for GDT, the unnormalised
    /// sum for TM).
    pub(crate) scores: Vec<f64>,
    /// The best transforms by the first slot, best first.
    pub(crate) leaders: Vec<Rigid>,
}

/// Searches for the superposition maximising `objective`, starting from the
/// known full-length fit `whole`.
pub(crate) fn maximise(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    objective: &Objective,
    whole: Rigid,
    schedule: Schedule,
) -> Outcome {
    let length = reference.len();
    let slots = objective.slots();
    let full = Problem { model, reference };
    let sampled = length > SAMPLE_CAP;

    let (sample_model, sample_reference);
    let search = if sampled {
        let picked: Vec<usize> = (0..SAMPLE_CAP).map(|j| j * length / SAMPLE_CAP).collect();
        sample_model = picked.iter().map(|&i| model[i]).collect::<Vec<_>>();
        sample_reference = picked.iter().map(|&i| reference[i]).collect::<Vec<_>>();
        Problem {
            model: &sample_model,
            reference: &sample_reference,
        }
    } else {
        full
    };

    let seed_list = seeds(search.reference.len());
    let rows = SeedRows::run(
        &seed_list,
        slots,
        schedule,
        || Scratch::new(search.reference.len(), slots),
        |scratch, seed, best| run_seed(search, objective, scratch, seed, best),
    );

    let mut total = Best::new(slots);
    if !sampled {
        let mut scratch = Scratch::new(length, slots);
        scratch.evaluate(full, objective, whole);
        total.offer(&scratch.scores, whole);
    }
    for seed in 0..rows.count() {
        total.view().merge(rows.scores(seed), rows.transforms(seed));
    }

    // The best seed transforms per slot, ties broken by seed order.
    let mut candidates: Vec<(usize, Rigid)> = Vec::new();
    for slot in 0..slots {
        let mut ranked: Vec<(f64, usize)> = (0..rows.count())
            .filter(|&seed| rows.transforms(seed)[slot].is_some())
            .map(|seed| (rows.scores(seed)[slot], seed))
            .collect();
        ranked.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        for &(_, index) in ranked.iter().take(CANDIDATES_PER_SLOT) {
            let held = rows.transforms(index)[slot];
            if let Some(transform) = held
                && !candidates.iter().any(|(_, known)| *known == transform)
            {
                candidates.push((slot, transform));
            }
        }
    }

    if !sampled {
        let mut leaders: Vec<(f64, Rigid)> = candidates
            .iter()
            .map(|&(slot, transform)| (rows.score_of(slot, transform), transform))
            .collect();
        leaders.push((total.scores[0], whole));
        return Outcome {
            scores: total.scores,
            leaders: ordered_leaders(leaders),
        };
    }

    polish_sampled(
        PolishInput {
            full,
            objective,
            whole,
            schedule,
        },
        &candidates,
        total,
    )
}

/// The inputs shared by every polish job of a sampled search.
#[derive(Clone, Copy)]
struct PolishInput<'a> {
    full: Problem<'a>,
    objective: &'a Objective,
    whole: Rigid,
    schedule: Schedule,
}

/// Takes every candidate (and the full-length fit, once per track) to the full
/// atom list, merging the polished results into `total`.
fn polish_sampled(
    input: PolishInput<'_>,
    candidates: &[(usize, Rigid)],
    mut total: Best,
) -> Outcome {
    let PolishInput {
        full,
        objective,
        whole,
        schedule,
    } = input;
    let length = full.reference.len();
    let slots = objective.slots();
    let mut jobs: Vec<(usize, Rigid)> = candidates
        .iter()
        .map(|&(slot, transform)| (slot, transform))
        .collect();
    for track in 0..objective.select.len() {
        jobs.push((track, whole));
    }
    let polish_job = |scratch: &mut Scratch, &(track_hint, transform): &(usize, Rigid)| {
        polish(
            full,
            objective,
            scratch,
            objective.select_for(track_hint),
            transform,
        )
    };
    let polished: Vec<Best> = match schedule {
        Schedule::Pool => jobs
            .par_iter()
            .map_init(|| Scratch::new(length, slots), polish_job)
            .collect(),
        Schedule::Inline => {
            let mut scratch = Scratch::new(length, slots);
            jobs.iter()
                .map(|job| polish_job(&mut scratch, job))
                .collect()
        }
    };
    let mut leaders = Vec::new();
    for result in &polished {
        total.merge(result);
        if let Some(transform) = result.transforms[0] {
            leaders.push((result.scores[0], transform));
        }
    }
    Outcome {
        scores: total.scores,
        leaders: ordered_leaders(leaders),
    }
}

/// Sorts leaders best first (stable), drops duplicates, and keeps a few.
fn ordered_leaders(mut leaders: Vec<(f64, Rigid)>) -> Vec<Rigid> {
    leaders.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mut out: Vec<Rigid> = Vec::new();
    for (_, transform) in leaders {
        if !out.contains(&transform) {
            out.push(transform);
        }
        if out.len() == CANDIDATES_PER_SLOT {
            break;
        }
    }
    out
}

/// Refines one transform on all atoms: score, keep the close atoms, refit.
fn polish(
    problem: Problem<'_>,
    objective: &Objective,
    scratch: &mut Scratch,
    cutoff: f64,
    start: Rigid,
) -> Best {
    let mut best = Best::new(objective.slots());
    // Score the transform itself first; it also yields the first fitting set.
    scratch.evaluate(problem, objective, start);
    best.offer(&scratch.scores, start);
    if !scratch.select_next(cutoff) {
        return best;
    }
    std::mem::swap(&mut scratch.current, &mut scratch.next);
    scratch.refine(
        problem,
        objective,
        cutoff,
        None,
        POLISH_ITERATIONS,
        &mut best.view(),
    );
    best
}

/// Normalises a count or sum by the atom count.
pub(crate) fn per_atom(value: f64, length: usize) -> f64 {
    value / usize_to_f64(length)
}

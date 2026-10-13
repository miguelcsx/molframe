//! What a seed search maximises, and the running best it keeps.

use molframe_geom::Rigid;

/// What is maximised.
#[derive(Clone, Debug)]
enum Kind {
    /// One count per cutoff: atoms within the cutoff.
    Counts(Vec<f64>),
    /// One sum of `1 / (1 + (d / d0)^2)`.
    Kernel(f64),
}

/// The quantity a search maximises and how its refinement picks atoms.
#[derive(Clone, Debug)]
pub(crate) struct Objective {
    kind: Kind,
    /// The distance under which an atom joins the next fitting set, per track.
    pub(crate) select: Vec<f64>,
}

impl Objective {
    /// The GDT objective: one tracked count and one refinement track per cutoff.
    pub(crate) fn gdt(cutoffs: &[f64]) -> Self {
        Self {
            kind: Kind::Counts(cutoffs.to_vec()),
            select: cutoffs.to_vec(),
        }
    }

    /// The TM objective with distance scale `d0`; one tracked sum, refined with
    /// two selection distances.
    pub(crate) fn tm(d0: f64) -> Self {
        Self {
            kind: Kind::Kernel(d0),
            select: vec![d0 + 1.5, d0.clamp(4.5, 8.0)],
        }
    }

    /// How many maxima are tracked.
    pub(crate) fn slots(&self) -> usize {
        match &self.kind {
            Kind::Counts(cutoffs) => cutoffs.len(),
            Kind::Kernel(_) => 1,
        }
    }

    pub(crate) fn score(&self, distances: &[f64], out: &mut [f64]) {
        match &self.kind {
            Kind::Counts(cutoffs) => {
                for (slot, cutoff) in out.iter_mut().zip(cutoffs) {
                    let mut count = 0u32;
                    for distance in distances {
                        count += u32::from(*distance <= *cutoff);
                    }
                    *slot = f64::from(count);
                }
            }
            Kind::Kernel(d0) => {
                let mut total = 0.0;
                for distance in distances {
                    let ratio = distance / d0;
                    total += 1.0 / (1.0 + ratio * ratio);
                }
                out[0] = total;
            }
        }
    }
}

/// The largest value seen for every tracked slot, with the transform that gave
/// it.
#[derive(Clone, Debug)]
pub(crate) struct Best {
    pub(crate) scores: Vec<f64>,
    pub(crate) transforms: Vec<Option<Rigid>>,
}

impl Best {
    pub(crate) fn new(slots: usize) -> Self {
        Self {
            scores: vec![-1.0; slots],
            transforms: vec![None; slots],
        }
    }

    pub(crate) fn offer(&mut self, scores: &[f64], transform: Rigid) {
        for ((best, held), &score) in self.scores.iter_mut().zip(&mut self.transforms).zip(scores) {
            if score > *best {
                *best = score;
                *held = Some(transform);
            }
        }
    }

    pub(crate) fn merge(&mut self, other: &Self) {
        for ((best, held), (&score, transform)) in self
            .scores
            .iter_mut()
            .zip(&mut self.transforms)
            .zip(other.scores.iter().zip(&other.transforms))
        {
            if score > *best {
                *best = score;
                *held = *transform;
            }
        }
    }
}

impl Objective {
    /// The selection cutoff for a candidate that came from slot `hint`.
    pub(crate) fn select_for(&self, hint: usize) -> f64 {
        self.select[hint.min(self.select.len() - 1)]
    }
}

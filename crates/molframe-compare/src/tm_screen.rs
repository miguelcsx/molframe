//! Opt-in screening of TM-score batches, so most pairs skip the seed search.
//!
//! At ten thousand structures an all-vs-all run is about 5e7 pairs, and almost
//! all of them are far below any threshold a caller cares about. The screen
//! drops those pairs cheaply and scores the rest exactly as [`crate::tm_score`]
//! does, so a surviving pair's value is bit-identical to the unfiltered one.
//!
//! Two stages run in order. The bound stage is exact: it prunes a pair only when
//! a proven upper bound on its TM-score is below the threshold, so no pair at or
//! above the threshold is ever lost. The optional coarse stage is a heuristic:
//! it prunes a pair whose single full-length fit scores below a fraction of the
//! threshold, which is fast but can lose a pair whose best superposition is
//! local (for example a multi-domain match), so its recall is a measured
//! property of the data, not a guarantee.
//!
//! Pruned pairs score `NaN`. Lengths are checked before filtering, so a length
//! mismatch is always an error; a degenerate fit is reported only for pairs
//! that reach scoring.

use crate::failure::CompareError;
use crate::seed_search::Schedule;
use crate::superposed::{tm_score_scheduled, tm_single_fit};
use crate::tm_batch::{condensed_pair, run};
use crate::tm_bound::TmBound;
use molframe_core::ExecutionContext;

/// How a screened pair was settled.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Stage {
    /// Scored by the full search without a screen.
    Scored,
    /// Scored by the full search after passing the screen.
    Evaluated,
    /// Dropped by the exact bound.
    Bound,
    /// Dropped by the coarse fit.
    Coarse,
}

/// Pair counts by outcome.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Tally {
    pub(crate) evaluated: u64,
    pub(crate) bound: u64,
    pub(crate) coarse: u64,
}

impl Tally {
    pub(crate) fn record(&mut self, stage: Stage) {
        match stage {
            Stage::Scored | Stage::Evaluated => self.evaluated += 1,
            Stage::Bound => self.bound += 1,
            Stage::Coarse => self.coarse += 1,
        }
    }

    pub(crate) fn add(&mut self, other: &Self) {
        self.evaluated += other.evaluated;
        self.bound += other.bound;
        self.coarse += other.coarse;
    }
}

/// The screen applied before scoring a pair.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TmPrefilter {
    threshold: f64,
    coarse_ratio: Option<f64>,
}

impl TmPrefilter {
    /// A screen that keeps every pair whose TM-score may reach `threshold`
    /// (exact, no coarse stage).
    ///
    /// # Errors
    ///
    /// Returns [`CompareError::InvalidScoreInput`] unless `threshold` is finite
    /// and in `(0, 1]`.
    pub fn new(threshold: f64) -> Result<Self, CompareError> {
        unit_interval(threshold)?;
        Ok(Self {
            threshold,
            coarse_ratio: None,
        })
    }

    /// Adds the heuristic coarse stage: drop a pair whose single-fit score is
    /// under `ratio * threshold`.
    ///
    /// # Errors
    ///
    /// Returns [`CompareError::InvalidScoreInput`] unless `ratio` is finite and
    /// in `(0, 1]`.
    pub fn with_coarse_ratio(mut self, ratio: f64) -> Result<Self, CompareError> {
        unit_interval(ratio)?;
        self.coarse_ratio = Some(ratio);
        Ok(self)
    }

    /// The score a pair must be able to reach to be kept.
    #[must_use]
    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    /// The coarse-stage ratio, when that stage is on.
    #[must_use]
    pub fn coarse_ratio(&self) -> Option<f64> {
        self.coarse_ratio
    }

    fn settle(
        &self,
        model: &[[f32; 3]],
        reference: &[[f32; 3]],
        bound: &mut TmBound,
    ) -> Result<(f64, Stage), CompareError> {
        if model.len() != reference.len() {
            return Err(CompareError::LengthMismatch {
                model: model.len(),
                reference: reference.len(),
            });
        }
        if model.len() >= 3 {
            if bound.bound(model, reference) < self.threshold {
                return Ok((f64::NAN, Stage::Bound));
            }
            if let Some(ratio) = self.coarse_ratio
                && let Ok(coarse) = tm_single_fit(model, reference)
                && coarse < ratio * self.threshold
            {
                return Ok((f64::NAN, Stage::Coarse));
            }
        }
        let score = tm_score_scheduled(model, reference, Schedule::Inline)?;
        Ok((score, Stage::Evaluated))
    }
}

fn unit_interval(value: f64) -> Result<(), CompareError> {
    if value.is_finite() && value > 0.0 && value <= 1.0 {
        Ok(())
    } else {
        Err(CompareError::InvalidScoreInput)
    }
}

/// Scores with `NaN` for pruned pairs, and how many pairs each stage settled.
#[derive(Clone, Debug)]
pub struct TmScreen {
    /// One value per requested pair; `NaN` when the screen dropped it.
    pub scores: Vec<f64>,
    /// Pairs that went through the full search.
    pub evaluated: u64,
    /// Pairs dropped by the exact bound.
    pub pruned_bound: u64,
    /// Pairs dropped by the coarse stage.
    pub pruned_coarse: u64,
}

/// [`crate::tm_score_pairs`] behind `prefilter`.
///
/// Every pair whose [`crate::tm_score`] is at least the threshold is scored
/// exactly; with the coarse stage off that is guaranteed, with it on it is
/// measured.
///
/// # Errors
///
/// As for [`crate::tm_score_pairs`].
pub fn tm_screen_pairs(
    structures: &[&[[f32; 3]]],
    pairs: &[(usize, usize)],
    prefilter: &TmPrefilter,
    context: &ExecutionContext,
) -> Result<TmScreen, CompareError> {
    if pairs
        .iter()
        .any(|&(a, b)| a >= structures.len() || b >= structures.len())
    {
        return Err(CompareError::InvalidMapping);
    }
    screen(structures, pairs.len(), prefilter, context, |index| {
        pairs.get(index).copied()
    })
}

/// [`crate::tm_score_all_pairs`] behind `prefilter`, in the same condensed
/// order.
///
/// # Errors
///
/// As for [`crate::tm_score_all_pairs`].
pub fn tm_screen_all_pairs(
    structures: &[&[[f32; 3]]],
    prefilter: &TmPrefilter,
    context: &ExecutionContext,
) -> Result<TmScreen, CompareError> {
    let count = structures.len();
    let total = count
        .checked_mul(count.saturating_sub(1))
        .map(|value| value / 2)
        .ok_or(CompareError::InvalidMapping)?;
    screen(structures, total, prefilter, context, |index| {
        condensed_pair(count, index)
    })
}

fn screen(
    structures: &[&[[f32; 3]]],
    total: usize,
    prefilter: &TmPrefilter,
    context: &ExecutionContext,
    pair_at: impl Fn(usize) -> Option<(usize, usize)> + Sync,
) -> Result<TmScreen, CompareError> {
    let longest = match structures.iter().map(|s| s.len()).max() {
        Some(longest) => longest,
        None => 0,
    };
    let (scores, tally) = run(
        structures,
        total,
        context,
        pair_at,
        TmBound::bytes(longest),
        |&(a, b), bound: &mut TmBound| prefilter.settle(structures[a], structures[b], bound),
    )?;
    Ok(TmScreen {
        scores,
        evaluated: tally.evaluated,
        pruned_bound: tally.bound,
        pruned_coarse: tally.coarse,
    })
}

#[cfg(test)]
#[path = "tm_screen_tests.rs"]
mod tests;

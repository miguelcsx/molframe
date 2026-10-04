//! Which decisions move the answer, alone and together.
//!
//! A plan is a full factorial: every combination of the varied decisions is
//! run once. That is enough to split the variation among the runs into what each
//! decision explains by itself, what each *pair* of decisions explains beyond
//! the sum of their separate effects, and what is left to higher-order
//! combinations.
//!
//! The split is the distance-based analysis of variance (Anderson's PERMANOVA
//! decomposition) over squared distances, so it needs only the metric the audit
//! already has and works for sets, numbers, vectors, rankings and networks alike.
//! For a scalar result with the absolute-error metric it is exactly the usual
//! sum-of-squares decomposition. For a distance that is not Euclidean (Jaccard
//! is not) a share can come out slightly negative; it is reported as computed
//! rather than clamped, because clamping would conceal that the metric and the
//! additive model disagree.
//!
//! Shares are fractions of the total variation, so they are comparable across
//! analyses with different units. A decision that matters only in combination
//! with another has a small main share and a large interaction share, which is
//! precisely the pattern a one-decision-at-a-time sweep cannot show.

use crate::numeric::usize_to_f64;
use molframe_core::contract::PolicyField;

/// What one decision explains on its own.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MainEffect {
    /// The decision.
    pub field: PolicyField,
    /// Mean distance between runs that differ in this decision only.
    pub mean_change: f64,
    /// How many pairs of runs that mean is over.
    pub comparisons: usize,
    /// Fraction of the total variation explained by this decision alone.
    pub share: f64,
}

/// What two decisions explain beyond their separate effects.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interaction {
    /// The first decision.
    pub first: PolicyField,
    /// The second decision.
    pub second: PolicyField,
    /// Fraction of the total variation that depends on the pair jointly.
    pub share: f64,
}

/// The decomposition of the variation among the runs of a plan.
#[derive(Clone, Debug, PartialEq)]
pub struct Decomposition {
    /// The sum of squared distances between runs divided by the run count; zero
    /// when every run gave the same answer.
    pub total_variation: f64,
    /// Mean distance over every pair of runs.
    pub mean_distance: f64,
    /// Largest distance between any two runs.
    pub max_distance: f64,
    /// One entry per varied decision, in plan order.
    pub main_effects: Vec<MainEffect>,
    /// One entry per pair of varied decisions, in plan order.
    pub interactions: Vec<Interaction>,
    /// What neither the main effects nor the pairs explain: combinations of
    /// three or more decisions.
    pub higher_order: f64,
}

/// Sums of squared distances within the groups of one grouping of the runs.
struct Groups {
    sums: Vec<f64>,
    sizes: Vec<usize>,
}

impl Groups {
    fn new(sizes: Vec<usize>) -> Self {
        Self {
            sums: vec![0.0; sizes.len()],
            sizes,
        }
    }

    /// The within-group sum of squares: each group's pairwise squared
    /// distances divided by its size.
    fn within(&self) -> f64 {
        self.sums
            .iter()
            .zip(&self.sizes)
            .map(|(sum, &size)| sum / usize_to_f64(size.max(1)))
            .sum()
    }
}

/// How many levels each axis has, read from the coordinates.
fn levels_of(coordinates: &[Vec<usize>], axes: usize) -> Vec<usize> {
    (0..axes)
        .map(|axis| {
            coordinates
                .iter()
                .filter_map(|point| point.get(axis))
                .max()
                .map_or(1, |&top| top + 1)
        })
        .collect()
}

/// One grouping per axis: the runs that share a level of it.
fn single_groups(levels: &[usize], coordinates: &[Vec<usize>]) -> Vec<Groups> {
    levels
        .iter()
        .enumerate()
        .map(|(axis, &count)| {
            Groups::new(
                (0..count)
                    .map(|level| {
                        coordinates
                            .iter()
                            .filter(|point| point.get(axis) == Some(&level))
                            .count()
                    })
                    .collect(),
            )
        })
        .collect()
}

/// One grouping per pair of axes: the runs that share a level of both.
fn pair_groups(levels: &[usize], coordinates: &[Vec<usize>]) -> Vec<((usize, usize), Groups)> {
    let mut pairs = Vec::new();
    for first in 0..levels.len() {
        for second in (first + 1)..levels.len() {
            let sizes = (0..levels[first] * levels[second])
                .map(|cell| {
                    let (a, b) = (cell / levels[second], cell % levels[second]);
                    coordinates
                        .iter()
                        .filter(|point| {
                            point.get(first) == Some(&a) && point.get(second) == Some(&b)
                        })
                        .count()
                })
                .collect();
            pairs.push(((first, second), Groups::new(sizes)));
        }
    }
    pairs
}

/// The sums every share is made from, gathered in one pass over the pairs of runs.
struct Sums {
    total: f64,
    distance: f64,
    maximum: f64,
    /// Per axis: the sum of distances between runs that differ in that axis only,
    /// and how many such pairs there are.
    one_axis: Vec<(f64, usize)>,
    single: Vec<Groups>,
    pairs: Vec<((usize, usize), Groups)>,
}

fn accumulate(
    coordinates: &[Vec<usize>],
    levels: &[usize],
    mut distance: impl FnMut(usize, usize) -> f64,
) -> Sums {
    let axes = levels.len();
    let mut sums = Sums {
        total: 0.0,
        distance: 0.0,
        maximum: 0.0,
        one_axis: vec![(0.0, 0); axes],
        single: single_groups(levels, coordinates),
        pairs: pair_groups(levels, coordinates),
    };
    for first in 0..coordinates.len() {
        for second in (first + 1)..coordinates.len() {
            let metric = distance(first, second);
            let square = metric * metric;
            sums.total += square;
            sums.distance += metric;
            sums.maximum = sums.maximum.max(metric);
            let (left, right) = (&coordinates[first], &coordinates[second]);
            let differing: Vec<usize> = (0..axes)
                .filter(|&axis| left[axis] != right[axis])
                .collect();
            if let [axis] = differing.as_slice() {
                sums.one_axis[*axis].0 += metric;
                sums.one_axis[*axis].1 += 1;
            }
            for (axis, groups) in sums.single.iter_mut().enumerate() {
                if left[axis] == right[axis] {
                    groups.sums[left[axis]] += square;
                }
            }
            for ((a, b), groups) in &mut sums.pairs {
                if left[*a] == right[*a] && left[*b] == right[*b] {
                    groups.sums[left[*a] * levels[*b] + left[*b]] += square;
                }
            }
        }
    }
    sums
}

/// Decomposes the variation among the runs of a full-factorial plan.
///
/// `fields` names the varied decisions; `coordinates[run]` gives each run's
/// level in each of them, as the plan lays them out. `distance(i, j)` is the
/// metric distance between runs `i` and `j`.
pub fn decompose(
    fields: &[PolicyField],
    coordinates: &[Vec<usize>],
    distance: impl FnMut(usize, usize) -> f64,
) -> Decomposition {
    let runs = coordinates.len();
    let levels = levels_of(coordinates, fields.len());
    let sums = accumulate(coordinates, &levels, distance);
    let count = usize_to_f64(runs.max(1));
    let total_variation = sums.total / count;
    let share = |explained: f64| {
        if total_variation > 0.0 {
            explained / total_variation
        } else {
            0.0
        }
    };
    let main: Vec<f64> = sums
        .single
        .iter()
        .map(|groups| total_variation - groups.within())
        .collect();
    let main_effects: Vec<MainEffect> = fields
        .iter()
        .zip(&main)
        .zip(&sums.one_axis)
        .map(
            |((&field, &explained), &(change, comparisons))| MainEffect {
                field,
                mean_change: if comparisons == 0 {
                    0.0
                } else {
                    change / usize_to_f64(comparisons)
                },
                comparisons,
                share: share(explained),
            },
        )
        .collect();
    let interactions: Vec<Interaction> = sums
        .pairs
        .iter()
        .map(|((a, b), groups)| Interaction {
            first: fields[*a],
            second: fields[*b],
            share: share(total_variation - groups.within() - main[*a] - main[*b]),
        })
        .collect();
    let explained: f64 = main_effects.iter().map(|effect| effect.share).sum::<f64>()
        + interactions.iter().map(|effect| effect.share).sum::<f64>();
    let pair_count = usize_to_f64(runs * runs.saturating_sub(1) / 2);
    Decomposition {
        total_variation,
        mean_distance: if pair_count == 0.0 {
            0.0
        } else {
            sums.distance / pair_count
        },
        max_distance: sums.maximum,
        main_effects,
        interactions,
        higher_order: if total_variation > 0.0 {
            1.0 - explained
        } else {
            0.0
        },
    }
}

#[cfg(test)]
#[path = "design_tests.rs"]
mod tests;

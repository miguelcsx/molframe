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

use crate::class::UncertaintyClass;
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

/// What a decision, or a class of decisions, explains once its overlaps with the
/// others are shared out fairly.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Attribution<K> {
    /// The decision or class.
    pub key: K,
    /// Its Shapley share of the total variation. The shares sum to one.
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
    /// Whether every combination of the varied decisions was run. Only then do the
    /// main-effect and interaction shares add up; the Shapley shares do not need it.
    pub balanced: bool,
    /// Each decision's Shapley share of the total variation, which is defined and sums
    /// to one however the runs are constrained. Empty when there are too many decisions
    /// to enumerate their coalitions, or no variation to share.
    pub shapley: Vec<Attribution<PolicyField>>,
    /// The same attribution with the four classes of uncertainty as the players.
    pub by_class: Vec<Attribution<UncertaintyClass>>,
}

/// The most decisions whose coalitions are enumerated for the Shapley shares.
const SHAPLEY_MAX_DECISIONS: usize = 14;

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
    /// The distance of every pair of runs, `first < second`, in row-major order.
    distances: Vec<f64>,
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
        distances: Vec::with_capacity(coordinates.len() * coordinates.len().saturating_sub(1) / 2),
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
            sums.distances.push(metric);
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
    let balanced = is_balanced(&levels, runs);
    let (shapley, by_class) = if total_variation > 0.0 && fields.len() <= SHAPLEY_MAX_DECISIONS {
        shapley_shares(fields, coordinates, &sums.distances, total_variation)
    } else {
        (Vec::new(), Vec::new())
    };
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
        balanced,
        shapley,
        by_class,
    }
}

/// Whether the runs are the whole product of the levels.
fn is_balanced(levels: &[usize], runs: usize) -> bool {
    levels
        .iter()
        .try_fold(1_usize, |product, &count| product.checked_mul(count))
        == Some(runs)
}

/// Shapley shares of the total variation, for each decision and for each class.
///
/// For a set `S` of decisions, `c(S)` is the fraction of the total variation removed by
/// holding the levels of `S` fixed: runs that agree on `S` are one group, and `c(S)` is one
/// minus the within-group sum of squares over the total. `c` is the closed Sobol index of
/// `S`, defined for any set of runs, so the attribution needs no balanced design. A
/// decision's Shapley share averages its marginal gain `c(S + i) - c(S)` over every order in
/// which the decisions could be fixed, and the shares sum to `c(all) = 1`.
fn shapley_shares(
    fields: &[PolicyField],
    coordinates: &[Vec<usize>],
    distances: &[f64],
    total_variation: f64,
) -> (
    Vec<Attribution<PolicyField>>,
    Vec<Attribution<UncertaintyClass>>,
) {
    let axes = fields.len();
    let runs = coordinates.len();
    let subsets = 1_usize << axes;
    // Size of the group each run falls in, for every set of decisions held fixed.
    let mut group_size = vec![vec![0_usize; runs]; subsets];
    for (mask, sizes) in group_size.iter_mut().enumerate() {
        let mut counts: std::collections::HashMap<Vec<usize>, usize> =
            std::collections::HashMap::new();
        let keys: Vec<Vec<usize>> = coordinates
            .iter()
            .map(|point| {
                (0..axes)
                    .filter(|axis| (mask >> axis) & 1 == 1)
                    .map(|axis| point[axis])
                    .collect()
            })
            .collect();
        for key in &keys {
            *counts.entry(key.clone()).or_insert(0) += 1;
        }
        for (run, key) in keys.iter().enumerate() {
            sizes[run] = match counts.get(key) {
                Some(&count) => count,
                None => 1,
            };
        }
    }
    // Within-group sum of squares for every set: each pair of runs counts in every set
    // whose levels the two runs share, weighted by the size of their group.
    let mut within = vec![0.0_f64; subsets];
    let mut offset = 0;
    for first in 0..runs {
        for second in (first + 1)..runs {
            let metric = distances[offset];
            offset += 1;
            let square = metric * metric;
            let same = (0..axes)
                .filter(|&axis| coordinates[first][axis] == coordinates[second][axis])
                .fold(0_usize, |mask, axis| mask | (1 << axis));
            // Every subset of `same`, the empty set included.
            let mut subset = same;
            loop {
                within[subset] += square / usize_to_f64(group_size[subset][first]);
                if subset == 0 {
                    break;
                }
                subset = (subset - 1) & same;
            }
        }
    }
    let total = total_variation;
    let closed: Vec<f64> = within
        .iter()
        .map(|sum| ((total - sum) / total).clamp(0.0, 1.0))
        .collect();
    let by_field = coalition_shares(axes, |mask| closed[mask])
        .into_iter()
        .zip(fields)
        .map(|(share, &key)| Attribution { key, share })
        .collect();
    let classes: Vec<UncertaintyClass> = UncertaintyClass::ALL
        .into_iter()
        .filter(|class| {
            fields
                .iter()
                .any(|&field| UncertaintyClass::of(field) == *class)
        })
        .collect();
    let mask_of = |chosen: usize| {
        fields
            .iter()
            .enumerate()
            .filter(|(_, field)| {
                classes.iter().enumerate().any(|(index, class)| {
                    (chosen >> index) & 1 == 1 && UncertaintyClass::of(**field) == *class
                })
            })
            .fold(0_usize, |mask, (axis, _)| mask | (1 << axis))
    };
    let by_class = coalition_shares(classes.len(), |chosen| closed[mask_of(chosen)])
        .into_iter()
        .zip(&classes)
        .map(|(share, &key)| Attribution { key, share })
        .collect();
    (by_field, by_class)
}

/// The Shapley value of each of `players` for a coalition value function.
fn coalition_shares(players: usize, value: impl Fn(usize) -> f64) -> Vec<f64> {
    let factorial: Vec<f64> = (0..=players)
        .scan(1.0_f64, |running, count| {
            if count > 0 {
                *running *= usize_to_f64(count);
            }
            Some(*running)
        })
        .collect();
    (0..players)
        .map(|player| {
            let mut share = 0.0;
            for coalition in 0..(1_usize << players) {
                if (coalition >> player) & 1 == 1 {
                    continue;
                }
                let size = (0..players)
                    .filter(|bit| (coalition >> bit) & 1 == 1)
                    .count();
                let weight = factorial[size] * factorial[players - size - 1] / factorial[players];
                share += weight * (value(coalition | (1 << player)) - value(coalition));
            }
            share
        })
        .collect()
}

#[cfg(test)]
#[path = "design_tests.rs"]
mod tests;

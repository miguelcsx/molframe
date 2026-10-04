use crate::metric::jaccard_distance;
use crate::{AuditPlan, AuditReport, AuditRun, DimensionSensitivity, SensitiveItem};
use std::collections::{BTreeMap, BTreeSet};

use crate::numeric::usize_to_f64;

/// Runs a planned audit and compares individually identifiable result items.
///
/// `analyse` owns the scientific calculation. `items` projects each result into
/// the identities whose presence should be audited. Both are caller-supplied so
/// this crate varies policies rather than reimplementing any analysis.
///
/// # Errors
///
/// Returns the first error produced by `analyse`; later policy points are not run.
pub fn audit<R, I, E, A, P>(
    plan: &AuditPlan,
    mut analyse: A,
    items: P,
) -> Result<AuditReport<R, I>, E>
where
    I: Ord + Clone,
    A: FnMut(&molframe_core::contract::AnalysisPolicy) -> Result<R, E>,
    P: Fn(&R) -> BTreeSet<I>,
{
    let mut runs = Vec::with_capacity(plan.cost());
    let mut sets = Vec::with_capacity(plan.cost());
    for policy in &plan.policies {
        let result = analyse(policy)?;
        sets.push(items(&result));
        runs.push(AuditRun {
            policy: policy.clone(),
            result,
        });
    }

    let frequencies = frequencies(&sets);
    let union_size = frequencies.len();
    let invariant = frequencies
        .values()
        .filter(|&&count| count == sets.len())
        .count();
    let stability = if union_size == 0 {
        1.0
    } else {
        usize_to_f64(invariant) / usize_to_f64(union_size)
    };
    let sensitive_items = frequencies
        .into_iter()
        .filter_map(|(item, count)| {
            (count != sets.len()).then(|| SensitiveItem {
                present_in: sets
                    .iter()
                    .enumerate()
                    .filter_map(|(index, set)| set.contains(&item).then_some(index))
                    .collect(),
                item,
            })
        })
        .collect();
    let decomposition =
        plan.decompose(|first, second| jaccard_distance(&sets[first], &sets[second]));
    let dimensions = plan
        .fields
        .iter()
        .enumerate()
        .zip(&decomposition.main_effects)
        .map(|((axis, &field), effect)| {
            dimension_report(
                field,
                axis,
                &plan.coordinates,
                &sets,
                effect.mean_change,
                effect.share,
            )
        })
        .collect();
    Ok(AuditReport {
        runs,
        stability,
        sensitive_items,
        dimensions,
        interactions: decomposition.interactions,
        higher_order: decomposition.higher_order,
        total_variation: decomposition.total_variation,
    })
}

fn frequencies<I: Ord + Clone>(sets: &[BTreeSet<I>]) -> BTreeMap<I, usize> {
    let mut counts = BTreeMap::new();
    for set in sets {
        for item in set {
            *counts.entry(item.clone()).or_insert(0) += 1;
        }
    }
    counts
}

fn dimension_report<I: Ord + Clone>(
    field: molframe_core::contract::PolicyField,
    axis: usize,
    coordinates: &[Vec<usize>],
    sets: &[BTreeSet<I>],
    mean_change: f64,
    main_effect_share: f64,
) -> DimensionSensitivity<I> {
    let mut changed = BTreeSet::new();
    for first in 0..sets.len() {
        for second in (first + 1)..sets.len() {
            if differ_only_on(axis, &coordinates[first], &coordinates[second]) {
                changed.extend(sets[first].symmetric_difference(&sets[second]).cloned());
            }
        }
    }
    DimensionSensitivity {
        field,
        sensitive_items: changed.into_iter().collect(),
        mean_change,
        main_effect_share,
    }
}

fn differ_only_on(axis: usize, first: &[usize], second: &[usize]) -> bool {
    first.get(axis) != second.get(axis)
        && first
            .iter()
            .zip(second)
            .enumerate()
            .all(|(index, (left, right))| index == axis || left == right)
}

#[cfg(test)]
#[path = "engine_tests.rs"]
mod tests;

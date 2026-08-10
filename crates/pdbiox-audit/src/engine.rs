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
    A: FnMut(&pdbiox_core::contract::AnalysisPolicy) -> Result<R, E>,
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
    let dimensions = plan
        .fields
        .iter()
        .enumerate()
        .map(|(axis, &field)| dimension_report(field, axis, &plan.coordinates, &sets))
        .collect();
    Ok(AuditReport {
        runs,
        stability,
        sensitive_items,
        dimensions,
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
    field: pdbiox_core::contract::PolicyField,
    axis: usize,
    coordinates: &[Vec<usize>],
    sets: &[BTreeSet<I>],
) -> DimensionSensitivity<I> {
    let mut changed = BTreeSet::new();
    let mut total_loss = 0.0;
    let mut comparisons = 0usize;
    for first in 0..sets.len() {
        for second in (first + 1)..sets.len() {
            if !differ_only_on(axis, &coordinates[first], &coordinates[second]) {
                continue;
            }
            changed.extend(sets[first].symmetric_difference(&sets[second]).cloned());
            total_loss += jaccard_loss(&sets[first], &sets[second]);
            comparisons += 1;
        }
    }
    DimensionSensitivity {
        field,
        sensitive_items: changed.into_iter().collect(),
        mean_change: if comparisons == 0 {
            0.0
        } else {
            total_loss / usize_to_f64(comparisons)
        },
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

fn jaccard_loss<I: Ord>(first: &BTreeSet<I>, second: &BTreeSet<I>) -> f64 {
    let union = first.union(second).count();
    if union == 0 {
        0.0
    } else {
        1.0 - usize_to_f64(first.intersection(second).count()) / usize_to_f64(union)
    }
}

#[cfg(test)]
#[path = "engine_tests.rs"]
mod tests;

//! Aggregating single-structure audits over a whole set of subjects.
//!
//! A per-structure audit reports how stable one analysis is under a policy
//! space. Running that over many structures and pooling the results is what a
//! population-scale sensitivity study needs: the mean and range of stability
//! across the set, and, for each varied policy dimension, how much it moves the
//! result on average. The scientific analysis and the item projection stay with
//! the caller, exactly as the single audit leaves them, so this only pools what
//! the audit already produced.
//!
//! Deterministic: results depend only on the subjects and the plan, and
//! dimensions are pooled by policy field so their order does not matter. Cost is
//! the single-audit cost times the number of subjects.

use std::collections::BTreeSet;

use molframe_core::contract::{AnalysisPolicy, PolicyField};

use crate::AuditPlan;
use crate::engine::audit;
use crate::numeric::usize_to_f64;

/// The pooled stability of one policy dimension across the subjects.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BatchDimension {
    /// The policy field varied.
    pub field: PolicyField,
    /// Its mean per-subject change, averaged over every subject.
    pub mean_change: f64,
}

/// The pooled result of auditing a set of subjects under one plan.
#[derive(Clone, Debug, PartialEq)]
pub struct BatchAudit {
    /// How many subjects were audited.
    pub subjects: usize,
    /// The mean stability across the subjects.
    pub mean_stability: f64,
    /// The least stable subject's stability.
    pub min_stability: f64,
    /// The most stable subject's stability.
    pub max_stability: f64,
    /// The per-dimension pooled sensitivity, ordered by first appearance.
    pub dimensions: Vec<BatchDimension>,
}

/// Audits every subject under `plan` and pools the stability statistics.
///
/// `analyse` runs the scientific calculation for one subject under one policy;
/// `items` projects a result into the identities whose presence is audited. An
/// empty subject set is perfectly stable, having nothing to destabilise.
///
/// # Errors
///
/// Returns the first error `analyse` produces; later subjects are not run.
pub fn audit_batch<S, R, I, E, A, P>(
    plan: &AuditPlan,
    subjects: &[S],
    analyse: A,
    items: P,
) -> Result<BatchAudit, E>
where
    I: Ord + Clone,
    A: Fn(&S, &AnalysisPolicy) -> Result<R, E>,
    P: Fn(&R) -> BTreeSet<I>,
{
    if subjects.is_empty() {
        return Ok(BatchAudit {
            subjects: 0,
            mean_stability: 1.0,
            min_stability: 1.0,
            max_stability: 1.0,
            dimensions: Vec::new(),
        });
    }

    let mut stabilities = Vec::with_capacity(subjects.len());
    let mut pooled: Vec<(PolicyField, f64)> = Vec::new();
    for subject in subjects {
        let report = audit(plan, |policy| analyse(subject, policy), &items)?;
        stabilities.push(report.stability);
        for dimension in &report.dimensions {
            match pooled
                .iter_mut()
                .find(|(field, _)| *field == dimension.field)
            {
                Some(entry) => entry.1 += dimension.mean_change,
                None => pooled.push((dimension.field, dimension.mean_change)),
            }
        }
    }

    let count = usize_to_f64(stabilities.len());
    let mean_stability = stabilities.iter().sum::<f64>() / count;
    let min_stability = stabilities.iter().copied().fold(f64::INFINITY, f64::min);
    let max_stability = stabilities
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let dimensions = pooled
        .into_iter()
        .map(|(field, sum)| BatchDimension {
            field,
            mean_change: sum / count,
        })
        .collect();

    Ok(BatchAudit {
        subjects: subjects.len(),
        mean_stability,
        min_stability,
        max_stability,
        dimensions,
    })
}

#[cfg(test)]
#[path = "batch_tests.rs"]
mod tests;

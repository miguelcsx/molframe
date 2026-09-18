//! Frozen measurement and decision conventions for external benchmarks.

use crate::{Comparison, MissingVerdict, Verdict, VerdictProfile, VerdictRule, VerdictStatus};
use std::collections::BTreeMap;

/// Atom subset used by one compatibility metric.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum ProfileAtomSet {
    /// C-alpha atoms over the complete designed scaffold.
    FullScaffoldCAlpha,
    /// Protein N, C-alpha, C and O atoms in the motif.
    MotifBackboneWithOxygen,
    /// Protein N, C-alpha and C atoms in catalytic residues.
    CatalyticBackbone,
    /// All non-hydrogen atoms in catalytic residues.
    CatalyticHeavyAtoms,
    /// Ligand atoms and protein backbone atoms used for the clash check.
    LigandAndBackbone,
}

/// Fit applied before a metric is measured.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum ProfileAlignment {
    /// Fit and measure the same atom correspondence.
    MeasuredAtoms,
    /// Fit catalytic N/C-alpha/C, then measure catalytic heavy atoms.
    CatalyticBackbone,
    /// No fit; measure a distance in the prediction's coordinate frame.
    PredictionFrame,
}

/// Frozen definition of one metric expected by a profile.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ProfileMetric {
    /// Stable key in the decomposed metric map.
    pub name: &'static str,
    /// Atom subset used by the metric.
    pub atoms: ProfileAtomSet,
    /// Alignment convention applied before measurement.
    pub alignment: ProfileAlignment,
}

/// How candidate-level verdicts become a scaffold-level result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CandidateAggregation {
    /// A scaffold passes when any predicted candidate passes every rule.
    Any,
}

/// A benchmark profile, including measurement conventions and aggregation.
#[derive(Clone, Debug, PartialEq)]
pub struct CompatibilityProfile {
    verdict: VerdictProfile,
    metrics: &'static [ProfileMetric],
    aggregation: CandidateAggregation,
}

impl CompatibilityProfile {
    /// Exact, versioned identifier.
    #[must_use]
    pub fn id(&self) -> &str {
        self.verdict.id()
    }

    /// Frozen metric definitions in decision order.
    #[must_use]
    pub const fn metrics(&self) -> &'static [ProfileMetric] {
        self.metrics
    }

    /// Candidate aggregation used by the external benchmark.
    #[must_use]
    pub const fn aggregation(&self) -> CandidateAggregation {
        self.aggregation
    }

    /// Applies the candidate-level thresholds.
    #[must_use]
    pub fn decide_candidate(&self, metrics: &BTreeMap<Box<str>, f64>) -> Verdict {
        self.verdict.decide(metrics)
    }

    /// Applies the benchmark's scaffold-level aggregation.
    #[must_use]
    pub fn decide_candidates<'a>(
        &self,
        candidates: impl IntoIterator<Item = &'a BTreeMap<Box<str>, f64>>,
    ) -> CompatibilityVerdict {
        let candidate_verdicts: Vec<_> = candidates
            .into_iter()
            .map(|metrics| self.decide_candidate(metrics))
            .collect();
        let status = aggregate_any(&candidate_verdicts);
        CompatibilityVerdict {
            profile: self.id().into(),
            status,
            candidates: candidate_verdicts,
        }
    }
}

/// Scaffold-level result retaining every candidate decision.
#[derive(Clone, Debug, PartialEq)]
pub struct CompatibilityVerdict {
    /// Exact profile identifier.
    pub profile: Box<str>,
    /// Aggregated decision.
    pub status: VerdictStatus,
    /// Candidate decisions in caller-provided order.
    pub candidates: Vec<Verdict>,
}

const MOTIFBENCH_METRICS: &[ProfileMetric] = &[
    ProfileMetric {
        name: "rmsd",
        atoms: ProfileAtomSet::FullScaffoldCAlpha,
        alignment: ProfileAlignment::MeasuredAtoms,
    },
    ProfileMetric {
        name: "motif_rmsd",
        atoms: ProfileAtomSet::MotifBackboneWithOxygen,
        alignment: ProfileAlignment::MeasuredAtoms,
    },
];

const AME_METRICS: &[ProfileMetric] = &[
    ProfileMetric {
        name: "catalytic_heavy_atom_rmsd",
        atoms: ProfileAtomSet::CatalyticHeavyAtoms,
        alignment: ProfileAlignment::CatalyticBackbone,
    },
    ProfileMetric {
        name: "ligand_backbone_min_distance",
        atoms: ProfileAtomSet::LigandAndBackbone,
        alignment: ProfileAlignment::PredictionFrame,
    },
];

/// `MotifBench` V1.0 compatibility profile.
///
/// The strict inequalities and failure-on-missing behavior intentionally match
/// the published evaluation pipeline. Candidate aggregation is `any` because a
/// generated scaffold succeeds when at least one fitted sequence succeeds.
#[must_use]
pub fn motifbench_1_0() -> CompatibilityProfile {
    CompatibilityProfile {
        verdict: VerdictProfile::new(
            "motifbench-1.0",
            [
                VerdictRule {
                    metric: "rmsd".into(),
                    comparison: Comparison::LessThan(2.0),
                },
                VerdictRule {
                    metric: "motif_rmsd".into(),
                    comparison: Comparison::LessThan(1.0),
                },
            ],
            MissingVerdict::Fail,
        ),
        metrics: MOTIFBENCH_METRICS,
        aggregation: CandidateAggregation::Any,
    }
}

/// Atomic Motif Enzyme heavy-atom compatibility profile.
///
/// A candidate needs catalytic heavy-atom RMSD below 1.5 angstrom after a
/// catalytic N/C-alpha/C fit, and no ligand-backbone pair closer than 1.5
/// angstrom. A scaffold succeeds when any predicted candidate succeeds.
#[must_use]
pub fn ame_heavy_atom_1_0() -> CompatibilityProfile {
    CompatibilityProfile {
        verdict: VerdictProfile::new(
            "ame-heavy-atom-1.0",
            [
                VerdictRule {
                    metric: "catalytic_heavy_atom_rmsd".into(),
                    comparison: Comparison::LessThan(1.5),
                },
                VerdictRule {
                    metric: "ligand_backbone_min_distance".into(),
                    comparison: Comparison::AtLeast(1.5),
                },
            ],
            MissingVerdict::Fail,
        ),
        metrics: AME_METRICS,
        aggregation: CandidateAggregation::Any,
    }
}

fn aggregate_any(candidates: &[Verdict]) -> VerdictStatus {
    if candidates
        .iter()
        .any(|candidate| candidate.status == VerdictStatus::Pass)
    {
        VerdictStatus::Pass
    } else if candidates.is_empty()
        || candidates
            .iter()
            .any(|candidate| candidate.status == VerdictStatus::Indeterminate)
    {
        VerdictStatus::Indeterminate
    } else {
        VerdictStatus::Fail
    }
}

#[cfg(test)]
#[path = "profile_tests.rs"]
mod tests;

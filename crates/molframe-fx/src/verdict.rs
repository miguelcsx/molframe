use std::collections::BTreeMap;

/// Numeric comparison made by a frozen profile rule.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum Comparison {
    /// Metric must be strictly less than the threshold.
    LessThan(f64),
    /// Metric must be at most the threshold.
    AtMost(f64),
    /// Metric must be strictly greater than the threshold.
    GreaterThan(f64),
    /// Metric must be at least the threshold.
    AtLeast(f64),
    /// Metric must lie in the inclusive range.
    Between(f64, f64),
}

impl Comparison {
    fn accepts(self, value: f64) -> bool {
        match self {
            Self::LessThan(limit) => value < limit,
            Self::AtMost(limit) => value <= limit,
            Self::GreaterThan(limit) => value > limit,
            Self::AtLeast(limit) => value >= limit,
            Self::Between(low, high) => value >= low && value <= high,
        }
    }
}

/// What a profile does when a required metric is absent.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissingVerdict {
    /// Return an indeterminate verdict.
    Indeterminate,
    /// Treat absence as a failed rule.
    Fail,
}

/// One named threshold rule.
#[derive(Clone, Debug, PartialEq)]
pub struct VerdictRule {
    /// Metric name in decomposed output.
    pub metric: Box<str>,
    /// Required comparison.
    pub comparison: Comparison,
}

/// Immutable, named and versioned mapping from metrics to a decision.
#[derive(Clone, Debug, PartialEq)]
pub struct VerdictProfile {
    id: Box<str>,
    rules: Vec<VerdictRule>,
    missing: MissingVerdict,
}

impl VerdictProfile {
    /// Creates a frozen profile value. The identifier should include its version.
    #[must_use]
    pub fn new(
        id: impl Into<Box<str>>,
        rules: impl IntoIterator<Item = VerdictRule>,
        missing: MissingVerdict,
    ) -> Self {
        Self {
            id: id.into(),
            rules: rules.into_iter().collect(),
            missing,
        }
    }

    /// Stable profile identifier including version.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Frozen rules in declaration order.
    #[must_use]
    pub fn rules(&self) -> &[VerdictRule] {
        &self.rules
    }

    /// Applies this profile to decomposed numeric measurements.
    #[must_use]
    pub fn decide(&self, metrics: &BTreeMap<Box<str>, f64>) -> Verdict {
        let mut outcomes = Vec::with_capacity(self.rules.len());
        let mut missing = false;
        let mut failed = false;
        for rule in &self.rules {
            let value = metrics.get(rule.metric.as_ref()).copied();
            let passed = value.map(|value| rule.comparison.accepts(value));
            missing |= passed.is_none();
            failed |=
                passed == Some(false) || (passed.is_none() && self.missing == MissingVerdict::Fail);
            outcomes.push(RuleOutcome {
                metric: rule.metric.clone(),
                value,
                passed,
            });
        }
        let status = if failed {
            VerdictStatus::Fail
        } else if missing {
            VerdictStatus::Indeterminate
        } else {
            VerdictStatus::Pass
        };
        Verdict {
            profile: self.id.clone(),
            status,
            outcomes,
        }
    }
}

/// Outcome of one profile rule.
#[derive(Clone, Debug, PartialEq)]
pub struct RuleOutcome {
    /// Metric name.
    pub metric: Box<str>,
    /// Observed value, if present.
    pub value: Option<f64>,
    /// `None` when the metric was absent.
    pub passed: Option<bool>,
}

/// Final decision state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VerdictStatus {
    /// Every rule passed.
    Pass,
    /// At least one rule failed.
    Fail,
    /// A required metric was absent under an indeterminate policy.
    Indeterminate,
}

/// Decision plus every decomposed rule outcome.
#[derive(Clone, Debug, PartialEq)]
pub struct Verdict {
    /// Exact profile identifier.
    pub profile: Box<str>,
    /// Decision state.
    pub status: VerdictStatus,
    /// Individual rule outcomes.
    pub outcomes: Vec<RuleOutcome>,
}

#[cfg(test)]
#[path = "verdict_tests.rs"]
mod tests;

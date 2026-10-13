//! Bounded aggregation of per-residue findings.

use molframe_core::diagnostic::{Code, Diagnostic};
use std::collections::BTreeMap;

/// Residues named in the sample of a summary finding.
const SAMPLE: usize = 8;
/// Components named individually in a summary finding.
const TOP_COMPONENTS: usize = 8;

/// Counts of one kind of per-residue finding.
///
/// A large model can raise the same finding for tens of thousands of residues.
/// The tally keeps exact totals and a bounded sample in file order, so the
/// report stays one diagnostic regardless of structure size.
#[derive(Default)]
pub(super) struct Tally {
    residues: usize,
    amount: usize,
    by_component: BTreeMap<Box<str>, usize>,
    sample: Vec<String>,
}

impl Tally {
    /// Records one residue, with the number of atoms involved when meaningful.
    pub(super) fn record(
        &mut self,
        component: &str,
        residue: impl std::fmt::Display,
        amount: usize,
    ) {
        self.residues += 1;
        self.amount += amount;
        match self.by_component.get_mut(component) {
            Some(count) => *count += 1,
            None => {
                self.by_component.insert(component.into(), 1);
            }
        }
        if self.sample.len() < SAMPLE {
            self.sample.push(sample_entry(component, &residue, amount));
        }
    }

    /// Absorbs a tally of later residues, keeping file order in the sample.
    pub(super) fn merge(&mut self, later: Self) {
        self.residues += later.residues;
        self.amount += later.amount;
        for (component, count) in later.by_component {
            *self.by_component.entry(component).or_insert(0) += count;
        }
        let room = SAMPLE.saturating_sub(self.sample.len());
        self.sample.extend(later.sample.into_iter().take(room));
    }

    /// The single summary finding, or `None` when nothing was recorded.
    pub(super) fn into_diagnostic(
        self,
        code: Code,
        amount_label: Option<&'static str>,
    ) -> Option<Diagnostic> {
        if self.residues == 0 {
            return None;
        }
        let mut ranked: Vec<(&Box<str>, &usize)> = self.by_component.iter().collect();
        ranked.sort_by(|left, right| right.1.cmp(left.1).then_with(|| left.0.cmp(right.0)));
        let mut by_component: Vec<String> = ranked
            .iter()
            .take(TOP_COMPONENTS)
            .map(|(component, count)| format!("{component}x{count}"))
            .collect();
        if ranked.len() > TOP_COMPONENTS {
            by_component.push(format!("+{} more", ranked.len() - TOP_COMPONENTS));
        }
        let mut diagnostic = Diagnostic::new(code)
            .with_context("residues", self.residues.to_string())
            .with_context("components", self.by_component.len().to_string());
        if let Some(label) = amount_label {
            diagnostic = diagnostic.with_context(label, self.amount.to_string());
        }
        Some(
            diagnostic
                .with_context("by component", by_component.join(", "))
                .with_context("sample", self.sample.join(", ")),
        )
    }
}

fn sample_entry(component: &str, residue: &dyn std::fmt::Display, amount: usize) -> String {
    if amount == 0 {
        format!("{component}@{residue}")
    } else {
        format!("{component}@{residue} ({amount} missing)")
    }
}

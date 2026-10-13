//! What [`Structure::with_distinct_labels`] changes, for a writer to report.
//!
//! The distinct view is silent by construction; a writer that projects it owes
//! its caller the list of identifiers it rewrote, so the rewrite never reads as
//! the depositor's own naming.

use super::Structure;

/// The normalised identifiers a distinct-label view rewrites.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct LabelRewrites {
    /// Each renamed chain as `(original label, new label)`.
    pub chains: Vec<(String, String)>,
    /// How many residues took a shifted normalised sequence number.
    pub residues: usize,
}

impl LabelRewrites {
    /// Whether the view changes nothing.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.chains.is_empty() && self.residues == 0
    }
}

impl Structure {
    /// The rewrites [`Structure::with_distinct_labels`] would make, in the
    /// first model's topology.
    #[must_use]
    pub fn label_rewrites(&self) -> LabelRewrites {
        let Some(distinct) = self.with_distinct_labels() else {
            return LabelRewrites::default();
        };
        let chains = self
            .data()
            .chains()
            .zip(distinct.data().chains())
            .filter_map(|(before, after)| {
                let (before, after) = (before.label()?, after.label()?);
                (before != after).then(|| (before.to_owned(), after.to_owned()))
            })
            .collect();
        let residues = self
            .data()
            .residues()
            .zip(distinct.data().residues())
            .filter(|(before, after)| before.label_seq_id() != after.label_seq_id())
            .count();
        LabelRewrites { chains, residues }
    }
}

#[cfg(test)]
#[path = "label_rewrites_tests.rs"]
mod tests;

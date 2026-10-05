//! Ordered composition preserves native topology and per-atom annotations.

use super::Structure;

impl Structure {
    /// Concatenates structures, preserving input row order and chemistry columns.
    ///
    /// # Errors
    /// Returns findings for incompatible frames, annotations or invalid topology.
    pub fn merge(structures: &[Self]) -> Result<Self, crate::Findings> {
        let sources: Vec<_> = structures
            .iter()
            .map(|value| value.engine().clone())
            .collect();
        molframe_core::Structure::merge(&sources)
            .map(Self::from)
            .map_err(crate::Findings::from)
    }
}

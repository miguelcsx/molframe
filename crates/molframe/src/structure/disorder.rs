//! Alternate-conformation resolution delegates to the shared storage engine.

use super::Structure;
use molframe_core::{
    AtomSelection,
    contract::{Analysis, AnalysisPolicy},
};

impl Structure {
    /// Selects alternate conformations while retaining coverage and decisions.
    ///
    /// Model and hydrogen restrictions are separate from alternate locations.
    #[must_use]
    pub fn resolve_altlocs(&self, policy: &AnalysisPolicy) -> Analysis<AtomSelection> {
        self.engine().resolve_altlocs(policy)
    }
}

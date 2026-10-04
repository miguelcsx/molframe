//! Altloc resolution for callers that cannot carry a refusal.

use super::Structure;
use crate::contract::AnalysisPolicy;
use crate::diagnostic::{Code, Diagnostic};
use crate::selection::AtomSelection;

impl Structure {
    /// The atoms an altloc policy selects, with the findings that came with them.
    ///
    /// For callers that cannot do anything with a refusal: it becomes the
    /// diagnostic that explains it. Code that can carry an indeterminate result
    /// should use [`Structure::resolve_altlocs`] and keep the outcome.
    ///
    /// # Errors
    ///
    /// Returns the first finding, or [`Code::E3001`], when the policy selects
    /// no consistent set of atoms.
    pub fn resolved_atoms(
        &self,
        policy: &AnalysisPolicy,
    ) -> Result<(AtomSelection, Vec<Diagnostic>), Diagnostic> {
        let resolution = self.resolve_altlocs(policy);
        let warnings = resolution.warnings.clone();
        resolution
            .into_result()
            .map(|atoms| (atoms, warnings.clone()))
            .map_err(|_| match warnings.into_iter().next() {
                Some(finding) => finding,
                None => Diagnostic::new(Code::E3001),
            })
    }
}

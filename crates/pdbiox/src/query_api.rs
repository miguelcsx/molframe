//! Structure-facing selection API over the query and spatial crates.

use pdbiox_core::contract::AnalysisPolicy;
use pdbiox_core::diagnostic::Diagnostic;
use pdbiox_core::structure::Structure;
use pdbiox_query::{Evaluation, Groups, Query};

/// Selection methods implemented by immutable structure snapshots.
pub trait QueryStructure {
    /// Evaluates a precompiled query with optional named groups.
    ///
    /// # Errors
    ///
    /// Returns diagnostics when required data is unavailable or a geometric
    /// parameter is invalid.
    fn select(
        &self,
        query: &Query,
        policy: &AnalysisPolicy,
        groups: &Groups,
    ) -> Result<Evaluation, Vec<Diagnostic>>;

    /// Compiles and evaluates textual syntax in one call.
    ///
    /// Use [`Query::compile`] directly when the same query will run repeatedly.
    ///
    /// # Errors
    ///
    /// Returns syntax, semantic or evaluation diagnostics.
    fn select_text(
        &self,
        source: &str,
        policy: &AnalysisPolicy,
    ) -> Result<Evaluation, Vec<Diagnostic>>;
}

impl QueryStructure for Structure {
    fn select(
        &self,
        query: &Query,
        policy: &AnalysisPolicy,
        groups: &Groups,
    ) -> Result<Evaluation, Vec<Diagnostic>> {
        #[cfg(feature = "spatial")]
        {
            let resolver = pdbiox_spatial::StructureSpatial::new(
                self,
                policy,
                pdbiox_spatial::SpatialBackend::Auto,
            )
            .map_err(|finding| vec![finding])?;
            query.evaluate(self, policy, groups, Some(&resolver))
        }
        #[cfg(not(feature = "spatial"))]
        {
            query.evaluate(self, policy, groups, None)
        }
    }

    fn select_text(
        &self,
        source: &str,
        policy: &AnalysisPolicy,
    ) -> Result<Evaluation, Vec<Diagnostic>> {
        let query = Query::compile(source)?;
        self.select(&query, policy, &Groups::new())
    }
}

#[cfg(test)]
#[path = "query_api_tests.rs"]
mod tests;

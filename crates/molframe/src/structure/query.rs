//! Structure-facing selection API over the query and spatial crates.

use molframe_core::ExecutionContext;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::diagnostic::Findings;
use molframe_core::structure::Structure;
use molframe_query::{Evaluation, Groups, Query};

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
        context: &ExecutionContext,
    ) -> Result<Evaluation, Findings>;

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
        context: &ExecutionContext,
    ) -> Result<Evaluation, Findings>;
}

impl QueryStructure for Structure {
    fn select(
        &self,
        query: &Query,
        policy: &AnalysisPolicy,
        groups: &Groups,
        context: &ExecutionContext,
    ) -> Result<Evaluation, Findings> {
        #[cfg(feature = "spatial")]
        {
            let resolver = molframe_spatial::StructureSpatial::new(
                self,
                policy,
                molframe_spatial::SpatialBackend::Auto,
                context,
            )
            .map_err(Findings::from)?;
            query
                .evaluate(self, policy, groups, Some(&resolver))
                .map_err(Findings::from)
        }
        #[cfg(not(feature = "spatial"))]
        {
            query
                .evaluate(self, policy, groups, None)
                .map_err(Findings::from)
        }
    }

    fn select_text(
        &self,
        source: &str,
        policy: &AnalysisPolicy,
        context: &ExecutionContext,
    ) -> Result<Evaluation, Findings> {
        let query = Query::compile(source)?;
        self.select(&query, policy, &Groups::new(), context)
    }
}

#[cfg(test)]
#[path = "query_tests.rs"]
mod tests;

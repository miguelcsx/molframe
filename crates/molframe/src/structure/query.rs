//! Structure-facing selection API over the query and spatial crates.

use molframe_core::ExecutionContext;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::diagnostic::Findings;
use molframe_core::structure::Structure as CoreStructure;
use molframe_query::{Evaluation, Groups, Query};

/// Selection methods implemented by immutable structure snapshots.
pub trait QueryStructure {
    /// Evaluates a precompiled query under the default execution context.
    ///
    /// Named groups are empty here; [`QueryStructure::select_with_options`]
    /// carries them, together with an explicit execution context, for callers
    /// that govern resources themselves.
    ///
    /// # Errors
    ///
    /// Returns diagnostics when required data is unavailable or a geometric
    /// parameter is invalid.
    fn select_query(&self, query: &Query, policy: &AnalysisPolicy) -> Result<Evaluation, Findings> {
        self.select_with_options(query, policy, &Groups::new(), &ExecutionContext::default())
    }

    /// Compiles and evaluates textual syntax in one call.
    ///
    /// Use [`Query::compile`] directly when the same query will run repeatedly.
    ///
    /// # Errors
    ///
    /// Returns syntax, semantic or evaluation diagnostics.
    fn select_text(&self, source: &str, policy: &AnalysisPolicy) -> Result<Evaluation, Findings> {
        let query = Query::compile(source)?;
        self.select_query(&query, policy)
    }

    /// Evaluates a precompiled query with explicit named groups and
    /// execution context.
    ///
    /// # Errors
    ///
    /// Returns the same diagnostics as [`QueryStructure::select_text`].
    fn select_with_options(
        &self,
        query: &Query,
        policy: &AnalysisPolicy,
        groups: &Groups,
        context: &ExecutionContext,
    ) -> Result<Evaluation, Findings>;
}

impl QueryStructure for CoreStructure {
    fn select_with_options(
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
            let _ = context;
            query
                .evaluate(self, policy, groups, None)
                .map_err(Findings::from)
        }
    }
}

#[cfg(test)]
#[path = "query_tests.rs"]
mod tests;

/// The same selection API on the facade's own handle, so a caller never has to
/// step through [`crate::Structure::engine`] to evaluate a query.
impl QueryStructure for crate::Structure {
    fn select_with_options(
        &self,
        query: &Query,
        policy: &AnalysisPolicy,
        groups: &Groups,
        context: &ExecutionContext,
    ) -> Result<Evaluation, Findings> {
        QueryStructure::select_with_options(self.engine(), query, policy, groups, context)
    }
}

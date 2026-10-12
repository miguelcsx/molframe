//! Structure-facing selection API over the query and spatial crates.

use molframe_core::ExecutionContext;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::diagnostic::{Code, Diagnostic, Findings};
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
        self.select_text_with_context(source, policy, &ExecutionContext::default())
    }

    /// Compiles text and evaluates spatial predicates with the caller's resource context.
    /// # Errors
    /// Returns query diagnostics, resource failures or cancellation diagnostics.
    fn select_text_with_context(
        &self,
        source: &str,
        policy: &AnalysisPolicy,
        context: &ExecutionContext,
    ) -> Result<Evaluation, Findings> {
        let query = Query::compile(source)?;
        self.select_with_options(&query, policy, &Groups::new(), context)
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
        context
            .cancellation()
            .check()
            .map_err(|_| Findings::from(Diagnostic::new(Code::E1904)))?;
        #[cfg(feature = "spatial")]
        let result = {
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
        };
        #[cfg(not(feature = "spatial"))]
        let result = {
            query
                .evaluate(self, policy, groups, None)
                .map_err(Findings::from)
        };
        context
            .cancellation()
            .check()
            .map_err(|_| Findings::from(Diagnostic::new(Code::E1904)))?;
        result
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

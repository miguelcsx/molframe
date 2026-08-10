//! Public query compilation and evaluation entry points.

use super::{Evaluation, Groups, Query};
use crate::builder::Builder;
use crate::plan::{LogicalPlan, PhysicalQuery};
use crate::spatial::SpatialResolver;
use pdbiox_core::contract::AnalysisPolicy;
use pdbiox_core::diagnostic::Diagnostic;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;

impl Query {
    /// Compiles textual selection syntax into a reusable typed plan.
    ///
    /// # Errors
    ///
    /// Returns registered diagnostics for lexical, syntax and type errors.
    pub fn compile(source: &str) -> Result<Self, Vec<Diagnostic>> {
        let parsed = crate::parser::parse(source)?;
        Ok(Self {
            expr: parsed.expr,
            warnings: parsed.warnings,
        })
    }

    /// Creates a reusable query from the typed builder surface.
    #[must_use]
    pub fn from_builder(builder: Builder) -> Self {
        Self {
            expr: builder.into_expr(),
            warnings: Vec::new(),
        }
    }

    /// Returns the storage-independent plan produced by either input surface.
    #[must_use]
    pub fn logical_plan(&self) -> LogicalPlan {
        LogicalPlan(self.expr.clone())
    }

    /// Binds this query to a structure dictionary and analysis policy.
    #[must_use]
    pub fn plan(&self, structure: &Structure, policy: &AnalysisPolicy) -> PhysicalQuery {
        crate::plan::bind_expr(&self.expr, structure, policy, self.warnings.clone())
    }

    /// Evaluates against the whole first-model topology.
    ///
    /// # Errors
    ///
    /// Returns diagnostics when required data or a spatial resolver is absent.
    pub fn evaluate(
        &self,
        structure: &Structure,
        policy: &AnalysisPolicy,
        groups: &Groups,
        spatial: Option<&dyn SpatialResolver>,
    ) -> Result<Evaluation, Vec<Diagnostic>> {
        self.evaluate_in(
            structure,
            &AtomSelection::All(structure.atom_count()),
            policy,
            groups,
            spatial,
        )
    }

    /// Evaluates inside a current view; `global` escapes to the whole topology.
    ///
    /// # Errors
    ///
    /// Returns the same diagnostics as [`Query::evaluate`].
    pub fn evaluate_in(
        &self,
        structure: &Structure,
        universe: &AtomSelection,
        policy: &AnalysisPolicy,
        groups: &Groups,
        spatial: Option<&dyn SpatialResolver>,
    ) -> Result<Evaluation, Vec<Diagnostic>> {
        self.plan(structure, policy)
            .evaluate_in(structure, universe, policy, groups, spatial)
    }
}

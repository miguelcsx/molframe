//! Public query compilation and evaluation entry points.

use super::{Evaluation, Groups, Query};
use crate::builder::Builder;
use crate::execution::QueryFingerprint;
use crate::plan::{LogicalPlan, PhysicalQuery};
use crate::spatial::SpatialResolver;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;
use std::ops::{BitAnd, BitOr, Not};

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
            source: source.into(),
        })
    }

    /// Creates a reusable query from the typed builder surface.
    #[must_use]
    pub fn from_builder(builder: Builder) -> Self {
        let (expr, source) = builder.into_parts();
        Self {
            expr,
            warnings: Vec::new(),
            source,
        }
    }

    /// Canonical parser input for this immutable query.
    #[must_use]
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Stable identity of the normalized typed query plan.
    #[must_use]
    pub fn fingerprint(&self) -> QueryFingerprint {
        QueryFingerprint::of(&self.expr)
    }

    /// Selects atoms within `radius` of `target`.
    #[must_use]
    pub fn within(radius: f32, target: Self) -> Self {
        let source = format!("within {radius} of ({})", target.source);
        Self {
            expr: crate::ast::Expr::Geometric(crate::ast::GeometricExpr::Within {
                radius,
                target: Box::new(target.expr),
            }),
            warnings: target.warnings,
            source: source.into(),
        }
    }

    /// Expands the selected atoms to their complete residues.
    #[must_use]
    pub fn by_residue(target: Self) -> Self {
        let source = format!("byres ({})", target.source);
        Self {
            expr: crate::ast::Expr::ByResidue(Box::new(target.expr)),
            warnings: target.warnings,
            source: source.into(),
        }
    }

    /// Selects complete residues within `radius` of `target`.
    #[must_use]
    pub fn residues_within(radius: f32, target: Self) -> Self {
        Self::by_residue(Self::within(radius, target))
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

impl From<Builder> for Query {
    fn from(builder: Builder) -> Self {
        Self::from_builder(builder)
    }
}

impl BitAnd for Query {
    type Output = Self;

    fn bitand(self, right: Self) -> Self::Output {
        let source = format!("({}) and ({})", self.source, right.source);
        Self {
            expr: crate::ast::Expr::And(Box::new(self.expr), Box::new(right.expr)),
            warnings: merged_warnings(self.warnings, right.warnings),
            source: source.into(),
        }
    }
}

impl BitOr for Query {
    type Output = Self;

    fn bitor(self, right: Self) -> Self::Output {
        let source = format!("({}) or ({})", self.source, right.source);
        Self {
            expr: crate::ast::Expr::Or(Box::new(self.expr), Box::new(right.expr)),
            warnings: merged_warnings(self.warnings, right.warnings),
            source: source.into(),
        }
    }
}

impl Not for Query {
    type Output = Self;

    fn not(self) -> Self::Output {
        let source = format!("not ({})", self.source);
        Self {
            expr: crate::ast::Expr::Not(Box::new(self.expr)),
            warnings: self.warnings,
            source: source.into(),
        }
    }
}

fn merged_warnings(mut left: Vec<Diagnostic>, right: Vec<Diagnostic>) -> Vec<Diagnostic> {
    left.extend(right);
    left.sort_by_key(Diagnostic::code);
    left
}

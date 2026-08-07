//! Compiled query evaluation over immutable structures.

use crate::ast::{Expr, GeometricExpr, SameKey};
use crate::builder::Builder;
use crate::parser;
use crate::plan::{LogicalPlan, PhysicalExpr, PhysicalQuery};
use crate::spatial::{GeometricRequest, SpatialResolver};
use pdbiox_core::contract::AnalysisPolicy;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;
use std::collections::BTreeMap;

/// Named selections supplied to `group <name>` expressions.
pub type Groups = BTreeMap<Box<str>, AtomSelection>;

/// A selected atom set and non-fatal findings produced while resolving it.
#[derive(Clone, Debug)]
pub struct Evaluation {
    /// The selected atom indices.
    pub selection: AtomSelection,
    /// Ambiguities and absent symbols worth reporting.
    pub warnings: Vec<Diagnostic>,
}

/// A parsed and reusable selection plan.
#[derive(Clone, Debug)]
pub struct Query {
    expr: Expr,
    warnings: Vec<Diagnostic>,
}

impl Query {
    /// Compiles textual selection syntax into a reusable typed plan.
    ///
    /// # Errors
    ///
    /// Returns registered diagnostics for lexical, syntax and type errors.
    pub fn compile(source: &str) -> Result<Self, Vec<Diagnostic>> {
        let parsed = parser::parse(source)?;
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
    ///
    /// The physical plan resolves globs once, folds constants and orders cheap
    /// predicates before expensive ones. It can be reused over dense frames.
    #[must_use]
    pub fn plan(&self, structure: &Structure, policy: &AnalysisPolicy) -> PhysicalQuery {
        crate::plan::bind(
            &self.logical_plan(),
            structure,
            policy,
            self.warnings.clone(),
        )
    }

    /// Evaluates against the whole first-model topology.
    ///
    /// # Errors
    ///
    /// Returns diagnostics when the query requests unavailable data or a
    /// geometric predicate without a spatial resolver.
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

impl PhysicalQuery {
    /// Executes this structure-bound plan over the complete topology.
    ///
    /// # Errors
    ///
    /// Returns diagnostics when data required by a predicate is unavailable.
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

    /// Executes inside a current view; `global` logical fallbacks can escape it.
    ///
    /// # Errors
    ///
    /// Returns the same diagnostics as [`PhysicalQuery::evaluate`].
    pub fn evaluate_in(
        &self,
        structure: &Structure,
        universe: &AtomSelection,
        policy: &AnalysisPolicy,
        groups: &Groups,
        spatial: Option<&dyn SpatialResolver>,
    ) -> Result<Evaluation, Vec<Diagnostic>> {
        if let Some(atom) = universe.iter().find(|atom| *atom >= structure.atom_count()) {
            return Err(vec![
                Diagnostic::new(Code::E6009).with_context("atom", atom.to_string()),
            ]);
        }
        let mut context = Context {
            structure,
            policy,
            groups,
            spatial,
            warnings: self.warnings.clone(),
            all: AtomSelection::All(structure.atom_count()),
        };
        let selection = context
            .eval_physical(&self.expr, universe)
            .map_err(|finding| vec![finding])?;
        context.warnings.sort_by_key(Diagnostic::code);
        Ok(Evaluation {
            selection,
            warnings: context.warnings,
        })
    }
}

struct Context<'a> {
    structure: &'a Structure,
    policy: &'a AnalysisPolicy,
    groups: &'a Groups,
    spatial: Option<&'a dyn SpatialResolver>,
    warnings: Vec<Diagnostic>,
    all: AtomSelection,
}

impl Context<'_> {
    fn eval_physical(
        &mut self,
        expr: &PhysicalExpr,
        universe: &AtomSelection,
    ) -> Result<AtomSelection, Diagnostic> {
        match expr {
            PhysicalExpr::All => Ok(universe.clone()),
            PhysicalExpr::None => Ok(AtomSelection::Empty),
            PhysicalExpr::And(left, right) => {
                let first = self.eval_physical(left, universe)?;
                if first.is_empty() {
                    return Ok(first);
                }
                Ok(first.intersect(&self.eval_physical(right, &first)?))
            }
            PhysicalExpr::Or(left, right) => Ok(self
                .eval_physical(left, universe)?
                .union(&self.eval_physical(right, universe)?)),
            PhysicalExpr::Not(target) => {
                Ok(universe.difference(&self.eval_physical(target, universe)?))
            }
            PhysicalExpr::ResolvedMembership { column, symbols } => {
                crate::predicate::membership_symbols(self.structure, universe, *column, symbols)
            }
            PhysicalExpr::Logical(logical) => self.eval(logical, universe),
        }
    }

    fn eval(&mut self, expr: &Expr, universe: &AtomSelection) -> Result<AtomSelection, Diagnostic> {
        match expr {
            Expr::All => Ok(universe.clone()),
            Expr::None => Ok(AtomSelection::Empty),
            Expr::And(left, right) => {
                let first = self.eval(left, universe)?;
                if first.is_empty() {
                    return Ok(first);
                }
                Ok(first.intersect(&self.eval(right, universe)?))
            }
            Expr::Or(left, right) => Ok(self
                .eval(left, universe)?
                .union(&self.eval(right, universe)?)),
            Expr::Not(target) => Ok(universe.difference(&self.eval(target, universe)?)),
            Expr::Global(target) => {
                let all = self.all.clone();
                self.eval(target, &all)
            }
            Expr::ByResidue(target) => {
                let selected = self.eval(target, universe)?;
                Ok(crate::expand::same_residue(
                    self.structure,
                    universe,
                    &selected,
                ))
            }
            Expr::Same { key, target } => {
                let selected = self.eval(target, universe)?;
                self.same(key, universe, &selected)
            }
            Expr::Bonded { depth, target } => {
                let selected = self.eval(target, universe)?;
                crate::connectivity::bonded(self.structure, universe, &selected, *depth)
            }
            Expr::Geometric(geometric) => self.geometric(geometric, universe),
            Expr::Comparison {
                column,
                operator,
                value,
                absolute,
            } => crate::predicate::comparison(
                self.structure,
                universe,
                *column,
                *operator,
                *value,
                *absolute,
                self.policy,
            ),
            Expr::Membership { column, values } => {
                let selected = crate::predicate::membership(
                    self.structure,
                    universe,
                    *column,
                    values,
                    self.policy,
                )?;
                if selected.is_empty() {
                    self.warnings.push(Diagnostic::new(Code::W4003));
                }
                Ok(selected)
            }
            Expr::Group(name) => {
                if let Some(group) = self.groups.get(name) {
                    Ok(universe.intersect(group))
                } else {
                    self.warnings
                        .push(Diagnostic::new(Code::W4003).with_context("group", name.to_string()));
                    Ok(AtomSelection::Empty)
                }
            }
            Expr::Atom {
                segment,
                residue,
                name,
            } => crate::predicate::atom_selector(
                self.structure,
                universe,
                segment,
                *residue,
                name,
                self.policy,
            ),
            Expr::Macro(macro_name) => {
                crate::macros::macro_selection(self.structure, universe, *macro_name)
            }
            Expr::Chirality(_) => Err(missing("CCD stereochemistry")),
            Expr::Smarts(_) => Err(missing("substructure matcher")),
        }
    }

    fn same(
        &self,
        key: &SameKey,
        universe: &AtomSelection,
        selected: &AtomSelection,
    ) -> Result<AtomSelection, Diagnostic> {
        match key {
            SameKey::Residue => Ok(crate::expand::same_residue(
                self.structure,
                universe,
                selected,
            )),
            SameKey::Chain => Ok(crate::expand::same_chain(
                self.structure,
                universe,
                selected,
            )),
            SameKey::Model => {
                if selected.is_empty() {
                    Ok(AtomSelection::Empty)
                } else {
                    Ok(universe.clone())
                }
            }
            SameKey::Entity => Ok(crate::expand::same_entity(
                self.structure,
                universe,
                selected,
            )),
            SameKey::Column(column) => crate::predicate::same_column(
                self.structure,
                universe,
                selected,
                *column,
                self.policy,
            ),
            SameKey::Fragment => {
                crate::connectivity::same_fragment(self.structure, universe, selected)
            }
            SameKey::Segment => crate::predicate::same_column(
                self.structure,
                universe,
                selected,
                crate::ast::Column::SegmentId,
                self.policy,
            ),
        }
    }

    fn geometric(
        &mut self,
        expr: &GeometricExpr,
        universe: &AtomSelection,
    ) -> Result<AtomSelection, Diagnostic> {
        let Some(spatial) = self.spatial else {
            return Err(missing("spatial resolver"));
        };
        let request = match expr {
            GeometricExpr::Within { radius, target } => {
                let target = self.eval(target, universe)?;
                return spatial.resolve(GeometricRequest::Within {
                    universe,
                    target: &target,
                    radius: *radius,
                    include_target: true,
                });
            }
            GeometricExpr::Beyond { radius, target } => {
                let target = self.eval(target, universe)?;
                return spatial.resolve(GeometricRequest::Beyond {
                    universe,
                    target: &target,
                    radius: *radius,
                });
            }
            GeometricExpr::Around { radius, target } => {
                let target = self.eval(target, universe)?;
                return spatial.resolve(GeometricRequest::Within {
                    universe,
                    target: &target,
                    radius: *radius,
                    include_target: false,
                });
            }
            GeometricExpr::SphereZone { radius, target } => {
                let target = self.eval(target, universe)?;
                return spatial.resolve(GeometricRequest::SphereZone {
                    universe,
                    target: &target,
                    radius: *radius,
                });
            }
            GeometricExpr::SphereLayer {
                inner,
                outer,
                target,
            } => {
                let target = self.eval(target, universe)?;
                return spatial.resolve(GeometricRequest::SphereLayer {
                    universe,
                    target: &target,
                    inner: *inner,
                    outer: *outer,
                });
            }
            GeometricExpr::IsoLayer {
                inner,
                outer,
                target,
            } => {
                let target = self.eval(target, universe)?;
                return spatial.resolve(GeometricRequest::IsoLayer {
                    universe,
                    target: &target,
                    inner: *inner,
                    outer: *outer,
                });
            }
            GeometricExpr::CylinderZone { .. } | GeometricExpr::CylinderLayer { .. } => {
                return self.geometric_cylinder(expr, universe, spatial);
            }
            GeometricExpr::Point { point, radius } => GeometricRequest::Point {
                universe,
                point: *point,
                radius: *radius,
            },
        };
        spatial.resolve(request)
    }

    fn geometric_cylinder(
        &mut self,
        expr: &GeometricExpr,
        universe: &AtomSelection,
        spatial: &dyn SpatialResolver,
    ) -> Result<AtomSelection, Diagnostic> {
        match expr {
            GeometricExpr::CylinderZone {
                radius,
                z_max,
                z_min,
                target,
            } => {
                let target = self.eval(target, universe)?;
                spatial.resolve(GeometricRequest::CylinderZone {
                    universe,
                    target: &target,
                    radius: *radius,
                    z_max: *z_max,
                    z_min: *z_min,
                })
            }
            GeometricExpr::CylinderLayer {
                inner,
                outer,
                z_max,
                z_min,
                target,
            } => {
                let target = self.eval(target, universe)?;
                spatial.resolve(GeometricRequest::CylinderLayer {
                    universe,
                    target: &target,
                    inner: *inner,
                    outer: *outer,
                    z_max: *z_max,
                    z_min: *z_min,
                })
            }
            _ => Err(Diagnostic::new(Code::E9001)),
        }
    }
}

fn missing(data: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E4003).with_context("required", data)
}

#[cfg(test)]
#[path = "eval_tests.rs"]
mod tests;

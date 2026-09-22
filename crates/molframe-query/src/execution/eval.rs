//! Compiled query evaluation over immutable structures.

use crate::ast::{Expr, GeometricExpr, SameKey};
use crate::plan::{PhysicalExpr, PhysicalQuery};
use crate::spatial::{GeometricRequest, SpatialResolver};
use molframe_chem::SmartsDataError;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;
use std::collections::BTreeMap;

#[path = "eval_query.rs"]
mod query;

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
    source: Box<str>,
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
        // The universe's own shape bounds every position it can name, so one
        // comparison replaces a walk over the atoms it contains. Asking about
        // each atom instead made every evaluation read the whole structure
        // before doing any work, which is what an `all` or `none` selection
        // cannot afford: their evaluation is otherwise constant.
        let atom_count = structure.atom_count();
        let bound = universe.position_bound();
        if bound > atom_count {
            return Err(vec![
                Diagnostic::new(Code::E6009).with_context("atom", bound.to_string()),
            ]);
        }

        let mut context = Context {
            structure,
            policy,
            groups,
            spatial,
            warnings: self.warnings.clone(),
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

/// Immutable evaluation dependencies plus evaluation-local warning state.
struct Context<'a> {
    structure: &'a Structure,
    policy: &'a AnalysisPolicy,
    groups: &'a Groups,
    spatial: Option<&'a dyn SpatialResolver>,
    warnings: Vec<Diagnostic>,
}

impl Context<'_> {
    /// Evaluates a lowered physical expression inside `universe`.
    ///
    /// Constant branches and pre-resolved memberships avoid falling back to the
    /// logical evaluator. Allocation is limited to selections produced by
    /// individual operators.
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

                let second = self.eval_physical(right, &first)?;
                Ok(first.intersect(&second))
            }
            PhysicalExpr::Or(left, right) => {
                let left = self.eval_physical(left, universe)?;
                let right = self.eval_physical(right, universe)?;

                Ok(left.union(&right))
            }
            PhysicalExpr::Not(target) => {
                let selected = self.eval_physical(target, universe)?;
                Ok(universe.difference(&selected))
            }
            PhysicalExpr::ResolvedMembership { column, symbols } => {
                crate::predicate::membership_symbols(self.structure, universe, *column, symbols)
            }
            PhysicalExpr::Logical(logical) => self.eval(logical, universe),
        }
    }

    /// Evaluates a logical expression inside `universe`.
    ///
    /// Boolean operations short-circuit empty conjunctions while preserving
    /// `global` semantics and the original diagnostic behavior.
    fn eval(&mut self, expr: &Expr, universe: &AtomSelection) -> Result<AtomSelection, Diagnostic> {
        match expr {
            Expr::All => Ok(universe.clone()),
            Expr::None => Ok(AtomSelection::Empty),
            Expr::And(left, right) => {
                let first = self.eval(left, universe)?;

                if first.is_empty() {
                    return Ok(first);
                }

                let second = self.eval(right, &first)?;
                Ok(first.intersect(&second))
            }
            Expr::Or(left, right) => {
                let left = self.eval(left, universe)?;
                let right = self.eval(right, universe)?;

                Ok(left.union(&right))
            }
            Expr::Not(target) => {
                let selected = self.eval(target, universe)?;
                Ok(universe.difference(&selected))
            }
            Expr::Global(target) => {
                let all = AtomSelection::All(self.structure.atom_count());
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
            Expr::Macro(macro_name) => crate::macros::macro_selection(
                self.structure,
                universe,
                *macro_name,
                &mut self.warnings,
            ),
            Expr::Chirality(configuration) => {
                crate::macros::chirality_selection(self.structure, universe, configuration)
            }
            Expr::Smarts(pattern) => smarts_selection(pattern, self.structure, universe),
        }
    }

    /// Expands `selected` according to a `same` relation.
    ///
    /// Delegates to specialized hierarchy, connectivity, or column expansion
    /// implementations while preserving `universe`.
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

    /// Resolves a geometric expression with the configured spatial backend.
    ///
    /// Target expressions are evaluated exactly once. Returns `E4003` when no
    /// spatial resolver is available.
    fn geometric(
        &mut self,
        expr: &GeometricExpr,
        universe: &AtomSelection,
    ) -> Result<AtomSelection, Diagnostic> {
        let Some(spatial) = self.spatial else {
            return Err(missing("spatial resolver"));
        };

        match expr {
            GeometricExpr::Within { radius, target } => {
                let target = self.eval(target, universe)?;

                spatial.resolve(GeometricRequest::Within {
                    universe,
                    target: &target,
                    radius: *radius,
                    include_target: true,
                })
            }
            GeometricExpr::Beyond { radius, target } => {
                let target = self.eval(target, universe)?;

                spatial.resolve(GeometricRequest::Beyond {
                    universe,
                    target: &target,
                    radius: *radius,
                })
            }
            GeometricExpr::Around { radius, target } => {
                let target = self.eval(target, universe)?;

                spatial.resolve(GeometricRequest::Within {
                    universe,
                    target: &target,
                    radius: *radius,
                    include_target: false,
                })
            }
            GeometricExpr::SphereZone { radius, target } => {
                let target = self.eval(target, universe)?;

                spatial.resolve(GeometricRequest::SphereZone {
                    universe,
                    target: &target,
                    radius: *radius,
                })
            }
            GeometricExpr::SphereLayer {
                inner,
                outer,
                target,
            } => {
                let target = self.eval(target, universe)?;

                spatial.resolve(GeometricRequest::SphereLayer {
                    universe,
                    target: &target,
                    inner: *inner,
                    outer: *outer,
                })
            }
            GeometricExpr::IsoLayer {
                inner,
                outer,
                target,
            } => {
                let target = self.eval(target, universe)?;

                spatial.resolve(GeometricRequest::IsoLayer {
                    universe,
                    target: &target,
                    inner: *inner,
                    outer: *outer,
                })
            }
            GeometricExpr::CylinderZone { .. } | GeometricExpr::CylinderLayer { .. } => {
                self.geometric_cylinder(expr, universe, spatial)
            }
            GeometricExpr::Point { point, radius } => spatial.resolve(GeometricRequest::Point {
                universe,
                point: *point,
                radius: *radius,
            }),
        }
    }

    /// Resolves cylinder-specific geometric expressions.
    ///
    /// This helper keeps the primary geometric dispatcher compact while
    /// evaluating the cylinder target exactly once.
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

/// Builds the standard missing-data diagnostic for an evaluation dependency.
///
/// Runtime and space are `O(1)` apart from diagnostic-owned context storage.
fn missing(data: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E4003).with_context("required", data)
}

fn missing_smarts_data(error: SmartsDataError) -> Diagnostic {
    let required = match error {
        SmartsDataError::Connectivity => "resolved connectivity",
        SmartsDataError::Aromaticity => "CCD aromaticity",
        SmartsDataError::FormalCharge => "CCD formal charges",
        SmartsDataError::Stereochemistry => "CCD stereochemistry",
    };
    missing(required)
}

fn smarts_selection(
    pattern: &molframe_chem::SmartsPattern,
    structure: &Structure,
    universe: &AtomSelection,
) -> Result<AtomSelection, Diagnostic> {
    let selected: AtomSelection = pattern
        .find_structure_matches(structure)
        .map_err(missing_smarts_data)?
        .into_iter()
        .flat_map(|matched| matched.atom_indices.into_vec())
        .filter_map(|atom| u32::try_from(atom).ok())
        .collect();
    Ok(universe.intersect(&selected))
}

#[cfg(test)]
#[path = "eval_tests.rs"]
mod tests;

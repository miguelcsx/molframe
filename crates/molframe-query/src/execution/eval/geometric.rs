//! Geometric predicate resolution for the evaluator.
//!
//! Evaluates each target expression exactly once, then hands the resolved
//! target to the configured spatial backend. Kept apart from the logical
//! evaluator because it is the only part that needs a spatial resolver.

use super::{Context, missing};
use crate::ast::GeometricExpr;
use crate::spatial::{GeometricRequest, SpatialResolver};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::selection::AtomSelection;

impl Context<'_> {
    /// Resolves a geometric expression with the configured spatial backend.
    ///
    /// Target expressions are evaluated exactly once. Returns `E4003` when no
    /// spatial resolver is available.
    pub(super) fn geometric(
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

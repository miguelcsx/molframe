//! Query-language request dispatch for a structure-bound spatial resolver.

use super::StructureSpatial;
use super::helpers::{resolve_beyond_request, resolve_within_request};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::selection::AtomSelection;
use molframe_query::{GeometricRequest, SpatialResolver};

impl SpatialResolver for StructureSpatial<'_> {
    /// Resolves one query-language geometric request.
    ///
    /// Each variant delegates to a focused helper to keep dispatch separate
    /// from individual spatial algorithms.
    fn resolve(&self, request: GeometricRequest<'_>) -> Result<AtomSelection, Diagnostic> {
        match request {
            GeometricRequest::Within {
                universe,
                target,
                radius,
                include_target,
            } => resolve_within_request(self, universe, target, radius, include_target),
            GeometricRequest::Beyond {
                universe,
                target,
                radius,
            } => resolve_beyond_request(self, universe, target, radius),
            GeometricRequest::SphereZone {
                universe,
                target,
                radius,
            } => self.radial_from_target(universe, target, 0.0, radius, None),
            GeometricRequest::SphereLayer {
                universe,
                target,
                inner,
                outer,
            } => self.radial_from_target(universe, target, inner, outer, None),
            GeometricRequest::IsoLayer {
                universe,
                target,
                inner,
                outer,
            } => self.iso_layer(universe, target, inner, outer),
            GeometricRequest::CylinderZone {
                universe,
                target,
                radius,
                z_max,
                z_min,
            } => self.radial_from_target(universe, target, 0.0, radius, Some((z_min, z_max))),
            GeometricRequest::CylinderLayer {
                universe,
                target,
                inner,
                outer,
                z_max,
                z_min,
            } => self.radial_from_target(universe, target, inner, outer, Some((z_min, z_max))),
            GeometricRequest::Point {
                universe,
                point,
                radius,
            } => self.radial(universe, point, 0.0, radius, None),
            _ => Err(Diagnostic::new(Code::E9001)),
        }
    }
}

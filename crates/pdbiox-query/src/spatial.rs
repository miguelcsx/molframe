//! Boundary between query planning and spatial execution.

use pdbiox_core::diagnostic::Diagnostic;
use pdbiox_core::selection::AtomSelection;

/// A fully planned geometric predicate.
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub enum GeometricRequest<'a> {
    /// Radius search around target atoms.
    Within {
        /// Atoms eligible for the result.
        universe: &'a AtomSelection,
        /// Centres of the radius search.
        target: &'a AtomSelection,
        /// Radius in ångström.
        radius: f32,
        /// Whether target atoms remain in the result.
        include_target: bool,
    },
    /// Atoms further than a radius from every target.
    Beyond {
        /// Atoms eligible for the result.
        universe: &'a AtomSelection,
        /// Centres excluded around.
        target: &'a AtomSelection,
        /// Exclusion radius in ångström.
        radius: f32,
    },
    /// Sphere around the target centre of geometry.
    SphereZone {
        /// Atoms eligible for the result.
        universe: &'a AtomSelection,
        /// Atoms defining the centre.
        target: &'a AtomSelection,
        /// Sphere radius in ångström.
        radius: f32,
    },
    /// Spherical shell around the target centre of geometry.
    SphereLayer {
        /// Atoms eligible for the result.
        universe: &'a AtomSelection,
        /// Atoms defining the centre.
        target: &'a AtomSelection,
        /// Inclusive inner radius.
        inner: f32,
        /// Inclusive outer radius.
        outer: f32,
    },
    /// Union of shells around each target atom.
    IsoLayer {
        /// Atoms eligible for the result.
        universe: &'a AtomSelection,
        /// Shell centres.
        target: &'a AtomSelection,
        /// Inclusive inner radius.
        inner: f32,
        /// Inclusive outer radius.
        outer: f32,
    },
    /// Cylinder around the target centre and global z axis.
    CylinderZone {
        /// Atoms eligible for the result.
        universe: &'a AtomSelection,
        /// Atoms defining the centre.
        target: &'a AtomSelection,
        /// Radial cutoff.
        radius: f32,
        /// Maximum signed z displacement.
        z_max: f32,
        /// Minimum signed z displacement.
        z_min: f32,
    },
    /// Cylindrical shell around the target centre and global z axis.
    CylinderLayer {
        /// Atoms eligible for the result.
        universe: &'a AtomSelection,
        /// Atoms defining the centre.
        target: &'a AtomSelection,
        /// Inclusive inner radial cutoff.
        inner: f32,
        /// Inclusive outer radial cutoff.
        outer: f32,
        /// Maximum signed z displacement.
        z_max: f32,
        /// Minimum signed z displacement.
        z_min: f32,
    },
    /// Sphere around an absolute Cartesian point.
    Point {
        /// Atoms eligible for the result.
        universe: &'a AtomSelection,
        /// Cartesian centre in ångström.
        point: [f32; 3],
        /// Sphere radius in ångström.
        radius: f32,
    },
}

/// Executes geometric plan nodes without making the query crate depend on a
/// particular spatial index implementation.
pub trait SpatialResolver {
    /// Evaluates one planned predicate.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when periodicity, coordinates or parameters make
    /// the operation undefined.
    fn resolve(&self, request: GeometricRequest<'_>) -> Result<AtomSelection, Diagnostic>;
}

//! Catalogued operations of the `surface` namespace.

use super::{Capability, Domain};

const D: Domain = Domain("surface");

pub(super) const ROWS: &[Capability] = &[
    D.row(
        "native_surface",
        "coordinates, radii, rows, provider",
        "MolecularSurface",
        "cpu,memory,temp",
    ),
    D.row("sasa", "coordinates, radii", "array", "cpu"),
    D.row("lee_richards", "coordinates, radii", "array", "cpu"),
    D.row("cavities", "coordinates, radii", "list", "cpu"),
    D.row(
        "surface_points",
        "coordinates, radii",
        "SurfacePoints",
        "cpu,memory",
    ),
    D.row(
        "surface_points_at_density",
        "coordinates, radii",
        "SurfacePoints",
        "cpu,memory",
    ),
    D.row(
        "solvent_excluded_surface",
        "coordinates, radii",
        "Mesh",
        "cpu,memory",
    ),
    D.row(
        "buried_surface",
        "coordinates, radii, in_first",
        "BuriedSurface",
        "cpu",
    ),
    D.row(
        "buried_solvent_excluded_surface",
        "coordinates, radii, roles",
        "BuriedSurface",
        "cpu,memory",
    ),
    D.row("atom_depths", "atoms, surface", "array", "cpu"),
];

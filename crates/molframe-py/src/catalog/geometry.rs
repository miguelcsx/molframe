//! Catalogued operations of the `geometry` namespace.

use super::{Capability, Domain};
use molframe::Cost;

const D: Domain = Domain("geometry");

pub(super) const ROWS: &[Capability] = &[
    D.row("centroid", "coordinates", "coordinate", "cpu")
        .cost(Cost::Borrow)
        .workflow(),
    D.row("distance_matrix", "coordinates", "matrix", "cpu,memory")
        .workflow(),
    D.row("rmsd", "coordinates, coordinates", "float", "cpu")
        .cost(Cost::Borrow)
        .workflow(),
    D.row("distances", "coordinates", "array", "cpu"),
    D.row("angles", "coordinates", "array", "cpu"),
    D.row("dihedrals", "coordinates", "array", "cpu"),
    D.row("centre_of_mass", "coordinates", "array", "cpu"),
    D.row("radius_of_gyration", "coordinates", "array", "cpu"),
    D.row("inertia_tensor", "coordinates", "array", "cpu"),
    D.row("principal_axes", "coordinates", "array", "cpu"),
    D.row("gyration_axes", "coordinates", "array", "cpu"),
    D.row("asphericity", "coordinates", "array", "cpu"),
    D.row("shape_parameter", "coordinates", "array", "cpu"),
    D.row("best_fit_plane", "coordinates", "array", "cpu"),
    D.row("plane_deviation", "coordinates", "array", "cpu"),
    D.row("rmsd_after_fit", "coordinates", "array", "cpu"),
    D.row("superpose", "coordinates", "array", "cpu"),
    D.row("rmsf", "coordinates", "array", "cpu"),
    D.row("backbone_torsions", "structure", "Table", "cpu"),
    D.row(
        "internal_coordinates",
        "structure",
        "InternalCoordinates",
        "cpu",
    ),
    D.row(
        "place_atom",
        "coordinates, internal coordinates",
        "coordinate",
        "cpu",
    ),
    D.row(
        "distance_matrix_blocks",
        "coordinates",
        "DistanceBlocks",
        "cpu,memory",
    ),
];

"""Numeric geometry kernels."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from .._native import geometry as _native

if TYPE_CHECKING:
    from collections.abc import Iterator

    from .. import ExecutionContext
    from . import Coordinates, DistanceMatrix

BatFrame = _native.BatFrame
InternalCoordinates = _native.InternalCoordinates
internal_coordinates = _native.internal_coordinates
place_atom = _native.place_atom
Superposition = _native.Superposition
centroid = _native.centroid
distance_matrix = _native.distance_matrix
rmsd = _native.rmsd
distances = _native.distances
angles = _native.angles
dihedrals = _native.dihedrals
centre_of_mass = _native.centre_of_mass
radius_of_gyration = _native.radius_of_gyration
inertia_tensor = _native.inertia_tensor
principal_axes = _native.principal_axes
gyration_axes = _native.gyration_axes
asphericity = _native.asphericity
shape_parameter = _native.shape_parameter
backbone_torsions = _native.backbone_torsions
best_fit_plane = _native.best_fit_plane
plane_deviation = _native.plane_deviation
rmsd_after_fit = _native.rmsd_after_fit
superpose = _native.superpose
rmsf = _native.rmsf


def distance_matrix_blocks(
    coordinates: Coordinates,
    block_rows: int | None = None,
    *,
    dtype: Literal["float64", "float32"] = "float64",
    context: ExecutionContext | None = None,
) -> Iterator[tuple[int, DistanceMatrix]]:
    """Yield ``(start_row, block)`` pairs covering the distance matrix.

    Each block has shape ``(rows, n)``; only one is resident at a time, so peak
    memory is bounded by ``block_rows * n`` values instead of ``n * n``. When
    ``block_rows`` is ``None`` it is chosen to fit the context memory budget.
    """
    size = coordinates.shape[0]
    if block_rows is None:
        block_rows = _native._distance_block_rows(  # noqa: SLF001  # pyright: ignore[reportPrivateUsage]
            coordinates, dtype=dtype, context=context
        )
    if block_rows < 1:
        message = "block_rows must be at least 1"
        raise ValueError(message)
    for start in range(0, size, block_rows):
        stop = min(start + block_rows, size)
        yield (
            start,
            _native._distance_rows(  # noqa: SLF001  # pyright: ignore[reportPrivateUsage]
                coordinates, start, stop, dtype=dtype, context=context
            ),
        )


__all__ = [
    "BatFrame",
    "InternalCoordinates",
    "Superposition",
    "angles",
    "asphericity",
    "backbone_torsions",
    "best_fit_plane",
    "centre_of_mass",
    "centroid",
    "dihedrals",
    "distance_matrix",
    "distance_matrix_blocks",
    "distances",
    "gyration_axes",
    "inertia_tensor",
    "internal_coordinates",
    "place_atom",
    "plane_deviation",
    "principal_axes",
    "radius_of_gyration",
    "rmsd",
    "rmsd_after_fit",
    "rmsf",
    "shape_parameter",
    "superpose",
]

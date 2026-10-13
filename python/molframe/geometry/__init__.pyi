from collections.abc import Sequence
from typing import Literal, Protocol

from .. import ExecutionContext, Float64Array, Structure, Table

class Coordinates(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class DistanceMatrix(Protocol):
    @property
    def shape(self) -> tuple[int, int]: ...

class Vector(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class Superposition:
    """A rigid transform ``R @ x + t`` and the deviation it achieves."""

    @property
    def rotation(self) -> Coordinates: ...
    @property
    def translation(self) -> list[float]: ...
    @property
    def rmsd(self) -> float: ...
    def apply(self, positions: Coordinates) -> Coordinates: ...

def centroid(coordinates: Coordinates) -> list[float] | None: ...
def distance_matrix(
    coordinates: Coordinates,
    *,
    dtype: Literal["float64", "float32"] = "float64",
    condensed: bool = False,
    context: ExecutionContext | None = None,
) -> DistanceMatrix:
    """All pairwise distances, charged to the context memory budget.

    ``dtype="float32"`` halves the memory; ``condensed=True`` returns the upper
    triangle as a 1-D array in ``scipy.spatial.distance.pdist`` order. Over
    budget it raises ``MemoryBudgetError`` naming these options.
    """

class DistanceBlocks:
    def __iter__(self) -> DistanceBlocks: ...
    def __next__(self) -> tuple[int, DistanceMatrix]: ...

def distance_matrix_blocks(
    coordinates: Coordinates,
    block_rows: int | None = None,
    *,
    dtype: Literal["float64", "float32"] = "float64",
    context: ExecutionContext | None = None,
) -> DistanceBlocks:
    """Stream ``(start_row, block)`` row blocks, one resident at a time.

    Each block has shape ``(rows, n)``, so peak memory is bounded by
    ``block_rows * n`` values instead of ``n * n``. When ``block_rows`` is
    ``None`` it is chosen to fit the context memory budget.
    """

def rmsd(mobile: Coordinates, reference: Coordinates) -> float: ...
def distances(left: Coordinates, right: Coordinates) -> Float64Array:
    """Distances between corresponding rows, in Å."""

def angles(first: Coordinates, vertex: Coordinates, third: Coordinates) -> Float64Array:
    """Angles at ``vertex`` in degrees; ``NaN`` where undefined."""

def dihedrals(
    first: Coordinates, second: Coordinates, third: Coordinates, fourth: Coordinates
) -> Float64Array:
    """Signed torsions in degrees in ``(-180, 180]``; ``NaN`` where undefined."""

def centre_of_mass(positions: Coordinates, masses: Vector | None = None) -> list[float] | None: ...
def radius_of_gyration(positions: Coordinates, masses: Vector | None = None) -> float | None: ...
def inertia_tensor(positions: Coordinates, masses: Vector | None = None) -> Coordinates | None: ...
def principal_axes(
    positions: Coordinates, masses: Vector | None = None
) -> tuple[Vector, Coordinates] | None:
    """Principal moments and axes (as columns), longest first."""

def gyration_axes(positions: Coordinates) -> tuple[Vector, Coordinates] | None: ...
def asphericity(positions: Coordinates) -> float | None: ...
def shape_parameter(positions: Coordinates) -> float | None: ...
def best_fit_plane(
    points: Coordinates,
) -> tuple[list[float], list[float]] | None:
    """Fit the plane through at least three points: ``(centre, normal)``."""

def plane_deviation(points: Coordinates) -> float | None: ...
def rmsd_after_fit(mobile: Coordinates, reference: Coordinates) -> float: ...
def superpose(mobile: Coordinates, reference: Coordinates) -> Superposition: ...
def rmsf(positions: Coordinates) -> Vector:
    """Per-atom fluctuation about the mean over ``(frames, atoms, 3)`` positions, in Å."""

def backbone_torsions(structure: Structure) -> Table:
    """Return phi, psi and omega in degrees for every protein residue (``NaN`` where undefined)."""

class BatFrame:
    """Native BAT frame, with angles and torsions in radians."""

    @property
    def seed_positions(self) -> list[list[float]]: ...
    @property
    def coordinates(self) -> list[list[float]]: ...

class InternalCoordinates:
    @property
    def seeds(self) -> list[tuple[int, list[float]]]: ...
    @property
    def atoms(
        self,
    ) -> list[tuple[int, list[int], list[float]]]:
        """Atom and reference indices; IJ length/angle/JK length/JKL angle/KL length/torsion."""

    def rebuild(self) -> list[list[float] | None]: ...
    def measure_bat(self, positions: Sequence[Sequence[float] | None]) -> BatFrame: ...
    def rebuild_bat(self, frame: BatFrame) -> list[list[float] | None]: ...

def internal_coordinates(structure: Structure, *, model: int = 0) -> InternalCoordinates:
    """Create a native deterministic bond-graph forest for a zero-based model."""

def place_atom(
    first: Sequence[float],
    second: Sequence[float],
    third: Sequence[float],
    length: float,
    angle: float,
    torsion: float,
) -> list[float] | None:
    """Place an atom using Å and radians; return None for degenerate references."""

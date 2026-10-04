from collections.abc import Sequence
from os import PathLike
from typing import Literal, Protocol

from .. import BoolArray, ExecutionContext, Float32Array, Float64Array, UInt32Array

class Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

def sasa(
    coordinates: Array,
    radii: Array,
    *,
    probe: float = 1.4,
    points: int = 960,
    context: ExecutionContext | None = None,
) -> Array: ...
def lee_richards(
    coordinates: Array,
    radii: Array,
    *,
    probe: float = 1.4,
    slices: int = 20,
    context: ExecutionContext | None = None,
) -> Array: ...
def cavities(
    coordinates: Array, radii: Array, *, probe: float = 1.4, resolution: float = 0.5
) -> list[tuple[float, tuple[float, float, float], int]]: ...

class UInt64Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class SurfacePoints:
    """Accessible-surface points: where, which way is out, and on which atom."""

    def __len__(self) -> int: ...
    @property
    def positions(self) -> Float32Array:
        """Point positions, ``(n, 3)``."""
    @property
    def normals(self) -> Float32Array:
        """Outward unit normals, ``(n, 3)``."""
    @property
    def atoms(self) -> UInt64Array:
        """The atom each point sits on."""

class BuriedSurface:
    """The accessible areas that decide an interface, and the area it buries."""

    @property
    def first_alone(self) -> float: ...
    @property
    def second_alone(self) -> float: ...
    @property
    def together(self) -> float: ...
    @property
    def buried(self) -> float:
        """``first_alone + second_alone - together``."""

class Curvatures:
    """Curvature descriptors at every vertex of a mesh."""

    def __len__(self) -> int: ...
    @property
    def mean(self) -> Float64Array: ...
    @property
    def gaussian(self) -> Float64Array: ...
    @property
    def maximum(self) -> Float64Array: ...
    @property
    def minimum(self) -> Float64Array: ...
    @property
    def shape_index(self) -> Float64Array: ...
    @property
    def curvedness(self) -> Float64Array: ...
    @property
    def quality(self) -> list[Literal["interior", "boundary", "degenerate"]]:
        """How far each estimate can be trusted."""

class Mesh:
    """An indexed triangle mesh, reported rather than repaired.

    It says what is wrong with it: boundary edges, edges shared by more than two faces and
    degenerate faces (including a face that names a vertex the mesh lacks). Distance and
    curvature operations refuse a mesh that is not manifold.
    """

    def __init__(self, vertices: Float32Array, faces: UInt32Array) -> None: ...
    @property
    def vertices(self) -> Float32Array: ...
    @property
    def faces(self) -> UInt32Array: ...
    @property
    def vertex_normals(self) -> Float32Array: ...
    @property
    def area(self) -> float:
        """Total area of the faces."""
    @property
    def is_manifold(self) -> bool: ...
    @property
    def boundary_edges(self) -> int: ...
    @property
    def non_manifold_edges(self) -> int: ...
    @property
    def degenerate_faces(self) -> int: ...
    @property
    def component_count(self) -> int:
        """Vertex-connected components."""
    def components(self) -> list[tuple[UInt32Array, float]]:
        """Face-connected components as ``(faces, area)``, in stable source-face order."""
    def filter_components(
        self, *, minimum_area: float = 0.0, maximum_components: int | None = None
    ) -> Mesh:
        """Keep components of at least ``minimum_area``, the largest ``maximum_components``."""
    def curvatures(self) -> Curvatures:
        """Curvature at every vertex; a non-manifold mesh is refused."""
    def geodesic_distances(self, source: int) -> Float64Array:
        """Distance from ``source`` along mesh edges (infinity marks another component)."""
    def patch(self, source: int, radius: float) -> UInt32Array:
        """Vertices within ``radius`` of ``source`` along mesh edges."""
    def write_obj(self, path: str | PathLike[str]) -> None:
        """Write Wavefront OBJ with vertex normals; an existing file is never replaced."""

def surface_points(
    coordinates: Array,
    radii: Array,
    *,
    probe: float = 1.4,
    samples: int = 960,
    context: ExecutionContext | None = None,
) -> SurfacePoints:
    """Return the accessible surface as points, ``samples`` test directions per atom."""

def surface_points_at_density(
    coordinates: Array,
    radii: Array,
    *,
    probe: float = 1.4,
    density: float = 4.0,
    context: ExecutionContext | None = None,
) -> SurfacePoints:
    """Return the accessible surface as points at ``density`` points per square angstrom."""

def solvent_excluded_surface(
    coordinates: Array,
    radii: Array,
    *,
    probe: float = 1.4,
    resolution: float = 0.5,
    max_cells: int | None = None,
    max_workspace_bytes: int | None = None,
) -> Mesh:
    """Return the solvent-excluded surface as a mesh, on a grid of cell edge ``resolution``."""

def buried_surface(
    coordinates: Array,
    radii: Array,
    in_first: BoolArray,
    *,
    probe: float = 1.4,
    points: int = 960,
    context: ExecutionContext | None = None,
) -> BuriedSurface:
    """Return the surface buried between the atoms marked ``in_first`` and all the others."""

def buried_solvent_excluded_surface(
    coordinates: Array,
    radii: Array,
    roles: Sequence[Literal["first", "second", "excluded"]],
    *,
    probe: float = 1.4,
    resolution: float = 0.5,
    max_cells: int | None = None,
    max_workspace_bytes: int | None = None,
) -> BuriedSurface:
    """Return the solvent-excluded surface buried between two molecules."""

def atom_depths(atoms: Array, surface: Array, *, cell_size: float = 3.0) -> Float32Array:
    """Return each atom's distance to the nearest surface point (infinity without any)."""

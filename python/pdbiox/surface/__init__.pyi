from pathlib import Path
from typing import final
from numpy import bool_, float32, float64
from numpy.typing import NDArray
from .._native.surface import *
from .._native.surface import __all__
from . import cavity, depth, geodesic, ses, slice_integration
from ._workflow import *
from ..core.contract import Analysis
from .._io_types import AtomDepthError, BuriedSurfaceError, SasaError, SurfaceGeometryError, SurfaceWorkflowError
from ..query import AnalysisPolicy
from .._trajectory import SurfaceMesh

class SurfacePoint:
    atom: int
    position: tuple[float, float, float]
    normal: tuple[float, float, float]

class ExcludedSurfacePoint:
    atom: int
    excluded: int
    position: tuple[float, float, float]
    area: float

class AtomContactArea:
    first: int
    second: int
    area: float

class BuriedSurface:
    first_alone: float
    second_alone: float
    together: float
    buried: float

@final
class AtomDepthOptions:
    def __init__(self, cell_size: float) -> None: ...

@final
class SurfaceGridOptions:
    STANDARD_MAX_CELLS: int
    def __init__(self, resolution: float, max_cells: int, max_workspace_bytes: int = ...) -> None: ...
    @staticmethod
    def standard(resolution: float) -> SurfaceGridOptions: ...

@final
class Cavity:
    volume: float
    representative: list[float]
    cells: int

class Sasa:
    def __init__(self, *, positions: NDArray[float32], radii: NDArray[float32], probe: float, sample_points: int) -> None: ...
    def execute(self) -> NDArray[float64]: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> Sasa: ...
    def __repr__(self) -> str: ...

class BuriedSurfaceOp:
    def __init__(self, *, positions: NDArray[float32], radii: NDArray[float32], first: NDArray[bool_], probe: float, sample_points: int) -> None: ...
    def execute(self) -> BuriedSurface: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> BuriedSurfaceOp: ...
    def __repr__(self) -> str: ...

class MoleculeRole:
    First: MoleculeRole
    Second: MoleculeRole
    Excluded: MoleculeRole

class SurfaceTriangle:
    vertices: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
    normal: tuple[float, float, float]
    area: float

class SolventExcludedSurface:
    triangles: list[SurfaceTriangle]
    area: float
    def indexed_mesh(self) -> IndexedSurfaceMesh: ...

class SurfaceFace:
    indices: tuple[int, int, int]
    def __init__(self, indices: tuple[int, int, int]) -> None: ...

class MeshReport:
    boundary_edges: int
    non_manifold_edges: int
    degenerate_faces: int
    components: int

class IndexedSurfaceMesh:
    def __init__(self, vertices: NDArray[float32], faces: list[tuple[int, int, int]]) -> None: ...
    vertices: NDArray[float32]
    faces: NDArray
    vertex_normals: NDArray[float32]
    report: MeshReport
    def neighbours(self, vertex: int) -> list[int]: ...

class SurfaceDistances:
    source: int
    distances: list[float]

class CurvatureQuality:
    Interior: CurvatureQuality
    Boundary: CurvatureQuality
    Degenerate: CurvatureQuality

class SurfaceCurvature:
    mean: float
    gaussian: float
    maximum: float
    minimum: float
    shape_index: float
    curvedness: float
    quality: CurvatureQuality

def fibonacci_sphere(count: int) -> list[tuple[float, float, float]]: ...
def shrake_rupley(positions: NDArray[float32], radii: NDArray[float32], probe: float, points: int) -> list[float]: ...
def solvent_accessible_surface(positions: NDArray[float32], radii: NDArray[float32], probe_radius: float, sample_points: int) -> list[float]: ...
def buried_surface(positions: NDArray[float32], radii: NDArray[float32], first: NDArray[bool_], probe_radius: float, sample_points: int) -> BuriedSurface: ...
def atom_depths(atoms: NDArray[float32], surface: NDArray[float32], options: AtomDepthOptions) -> NDArray[float32]: ...
def cavities(positions: NDArray[float32], radii: NDArray[float32], probe: float, options: SurfaceGridOptions) -> list[Cavity]: ...
def lee_richards(positions: NDArray[float32], radii: NDArray[float32], probe: float, slices: int) -> list[float]: ...
def surface_points(positions: NDArray[float32], radii: NDArray[float32], probe: float, samples: int) -> list[SurfacePoint]: ...
def surface_points_at_density(positions: NDArray[float32], radii: NDArray[float32], probe: float, density: float) -> list[SurfacePoint]: ...
def atom_contact_areas(positions: NDArray[float32], radii: NDArray[float32], probe: float, density: float) -> list[AtomContactArea]: ...
def surface_points_excluding_pairs(positions: NDArray[float32], radii: NDArray[float32], probe: float, density: float, pairs: list[tuple[int, int]]) -> list[ExcludedSurfacePoint]: ...
def buried_solvent_excluded_surface(positions: NDArray[float32], radii: NDArray[float32], roles: list[MoleculeRole], probe: float, resolution: float) -> BuriedSurface: ...
def buried_solvent_excluded_surface_with_options(positions: NDArray[float32], radii: NDArray[float32], roles: list[MoleculeRole], probe: float, options: SurfaceGridOptions) -> BuriedSurface: ...
def cavities_with_options(positions: NDArray[float32], radii: NDArray[float32], probe: float, options: SurfaceGridOptions) -> list[Cavity]: ...
def solvent_excluded_surface(positions: NDArray[float32], radii: NDArray[float32], probe: float, resolution: float) -> SolventExcludedSurface: ...
def solvent_excluded_surface_with_options(positions: NDArray[float32], radii: NDArray[float32], probe: float, options: SurfaceGridOptions) -> SolventExcludedSurface: ...
def edge_geodesic_distances(mesh: IndexedSurfaceMesh, source: int) -> SurfaceDistances: ...
def surface_patch(mesh: IndexedSurfaceMesh, source: int, radius: float) -> list[int]: ...
def surface_curvatures(mesh: IndexedSurfaceMesh) -> list[SurfaceCurvature]: ...
def write_obj(path: str | Path, mesh: IndexedSurfaceMesh) -> None: ...
def governed_surface_geometry(positions: NDArray[float32], radii: NDArray[float32], options: SurfaceWorkflowOptions, policy: AnalysisPolicy) -> Analysis[SurfaceWorkflowResult]: ...

@final
class SurfaceComponent:
    faces: list[int]
    area: float

@final
class SurfaceComponentFilter:
    def __init__(self, minimum_area: float = ..., maximum_components: int | None = ...) -> None: ...
    @property
    def minimum_area(self) -> float: ...
    @property
    def maximum_components(self) -> int | None: ...

@final
class SurfaceComponentError:
    InvalidFilter: SurfaceComponentError
    MeshTooLarge: SurfaceComponentError

def surface_components(mesh: IndexedSurfaceMesh) -> list[SurfaceComponent]: ...
def filter_surface_components(mesh: IndexedSurfaceMesh, filter: SurfaceComponentFilter) -> IndexedSurfaceMesh: ...

from collections.abc import Callable
from ..core.execution import ExecutionContext
from ..xtal import UnitCell

def visit_shrake_rupley(positions: NDArray[float32], radii: NDArray[float32], probe: float, samples: int, emit: Callable[[int, float], None], context: ExecutionContext, *, cell: UnitCell | None = None) -> None: ...

def collect_shrake_rupley(positions: NDArray[float32], radii: NDArray[float32], probe: float, samples: int, context: ExecutionContext, *, cell: UnitCell | None = None) -> NDArray[float64]: ...

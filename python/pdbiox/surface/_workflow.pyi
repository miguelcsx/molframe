from typing import final
from numpy import float32, float64
from numpy.typing import NDArray
from ..core.contract import Analysis
from ..query import AnalysisPolicy
from .._trajectory import SurfaceMesh
from . import SurfaceGridOptions

@final
class SurfaceWorkflowResult:
    mesh: SurfaceMesh
    curvature: list[tuple[float, float, float, float, float, float, str]]
    geodesic_source: int | None
    geodesic_distances: NDArray[float64] | None

@final
class SurfaceWorkflowOptions:
    def __init__(self, selection: str, radii_set: str, probe: float, grid: SurfaceGridOptions, source_vertex: int | None = None) -> None: ...

def analyse_surface_geometry(positions: NDArray[float32], radii: NDArray[float32], options: SurfaceWorkflowOptions, policy: AnalysisPolicy) -> Analysis[SurfaceWorkflowResult]: ...

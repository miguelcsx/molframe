from os import PathLike
from typing import final
from numpy import float32, float64, int64, uint32
from numpy.typing import NDArray
from numpy.ma import MaskedArray
from .surface import SurfaceGridOptions
@final
class CartesianFit:
    @staticmethod
    def none() -> CartesianFit: ...
    @staticmethod
    def reference(coordinates: NDArray[float32]) -> CartesianFit: ...
    @staticmethod
    def iterative_mean(max_iterations: int, tolerance: float) -> CartesianFit: ...
from ._trajectory_dispatch import AmberAsciiReadOptions, TrajectoryReadOptions, read_trajectory, write_trajectory
DEFAULT_PAIRWISE_MEMORY_LIMIT: int
DEFAULT_PROCRUSTES_TOLERANCE: float
DEFAULT_PROCRUSTES_MAXIMUM_ITERATIONS: int

from ._trajectory_format_models import *
from ._trajectory_runtime import *
from ._trajectory_operations import *

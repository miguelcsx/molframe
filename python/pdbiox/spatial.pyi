from typing import final
from numpy import float32
from numpy.typing import NDArray
from .crystallography import UnitCell
from .graph import SpatialBackend

@final
class NeighborPair:
    first: int
    second: int
    distance: float

@final
class AutoBackendProfile:
    def __init__(self, brute_force_pair_limit: int, kd_target_minimum: int, kd_query_ratio: int, periodic_backend: SpatialBackend) -> None: ...
    @staticmethod
    def balanced() -> AutoBackendProfile: ...
    brute_force_pair_limit: int
    kd_target_minimum: int
    kd_query_ratio: int
    periodic_backend: SpatialBackend

@final
class NeighborSkinProfile:
    def __init__(self, cutoff_ratio: float, minimum: float) -> None: ...
    @staticmethod
    def balanced() -> NeighborSkinProfile: ...
    cutoff_ratio: float
    minimum: float

@final
class CellGridOptions:
    def __init__(self, maximum_cell_count: int, edge_growth_factor: float) -> None: ...
    @staticmethod
    def memory_balanced() -> CellGridOptions: ...
    maximum_cell_count: int
    edge_growth_factor: float

@final
class KdPeriodicOptions:
    def __init__(self, maximum_image_count: int) -> None: ...
    @staticmethod
    def balanced() -> KdPeriodicOptions: ...
    maximum_image_count: int

@final
class SpatialSearchOptions:
    def __init__(self, backend: SpatialBackend = SpatialBackend.Auto, *, automatic: AutoBackendProfile | None = None, neighbor_skin: NeighborSkinProfile | None = None, cell_grid: CellGridOptions | None = None, kd_periodic: KdPeriodicOptions | None = None) -> None: ...
    @staticmethod
    def balanced() -> SpatialSearchOptions: ...
    def plan(self, left_count: int, right_count: int, periodic: bool, cutoff: float) -> SpatialPlan: ...

@final
class SpatialPlan:
    requested_backend: SpatialBackend
    backend: SpatialBackend
    neighbor_skin: float | None
    options: SpatialSearchOptions

def neighbor_pairs(positions: NDArray[float32], cutoff: float, *, left: list[int] | None = None, right: list[int] | None = None, backend: SpatialBackend = SpatialBackend.Auto, cell: UnitCell | None = None) -> list[NeighborPair]: ...
def atoms_within(positions: NDArray[float32], targets: list[int], cutoff: float, *, query: list[int] | None = None, backend: SpatialBackend = SpatialBackend.Auto, cell: UnitCell | None = None) -> list[int]: ...
def neighbor_pairs_with_options(positions: NDArray[float32], cutoff: float, options: SpatialSearchOptions, *, left: list[int] | None = None, right: list[int] | None = None, cell: UnitCell | None = None) -> list[NeighborPair]: ...
def atoms_within_with_options(positions: NDArray[float32], targets: list[int], cutoff: float, options: SpatialSearchOptions, *, query: list[int] | None = None, cell: UnitCell | None = None) -> list[int]: ...

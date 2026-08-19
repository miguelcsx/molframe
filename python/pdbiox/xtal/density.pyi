from typing import final
from numpy import float32
from numpy.typing import NDArray
from ..xtal import UnitCell
from .._io_types import MapStatisticsError, MrcError

@final
class MapBoundary:
    Missing: MapBoundary
    Periodic: MapBoundary

@final
class MapStatistics:
    count: int
    minimum: float
    maximum: float
    mean: float
    sigma: float

@final
class MapHistogram:
    minimum: float
    maximum: float
    counts: list[int]

@final
class DensityMap:
    def __init__(self, dimensions: tuple[int, int, int], starts: tuple[int, int, int], sampling: tuple[int, int, int], cell: UnitCell, origin: tuple[float, float, float], space_group: int, labels: list[str], extended_header: bytes, values: list[float]) -> None: ...
    @staticmethod
    def from_mrc_bytes(data: bytes) -> DensityMap: ...
    dimensions: tuple[int, int, int]
    starts: tuple[int, int, int]
    sampling: tuple[int, int, int]
    cell: UnitCell
    origin: tuple[float, float, float]
    space_group: int
    labels: list[str]
    extended_header: bytes
    values: NDArray[float32]
    def to_mrc_bytes(self) -> bytes: ...
    def statistics(self) -> MapStatistics: ...
    def masked_statistics(self, mask: list[bool]) -> MapStatistics: ...
    def histogram(self, bins: int, minimum: float, maximum: float) -> MapHistogram: ...

from typing import final
from numpy import float32
from numpy.typing import NDArray
from ..xtal import UnitCell
from .._io_types import MapStatisticsError, MrcError
from .._io_types import MrcBrickError
from .._provider import ChunkId, DatasetId

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

@final
class CubeAtom:
    number: int
    charge: float
    position: tuple[float, float, float]

@final
class CubeGrid:
    map: DensityMap
    atoms: list[CubeAtom]
    orbitals: list[int]
    fields: int

@final
class MrcBlockOptions:
    def __init__(self, *, memory_limit_bytes: int = ...) -> None: ...
    memory_limit_bytes: int

@final
class MrcMapDescriptor:
    dimensions: tuple[int, int, int]
    starts: tuple[int, int, int]
    sampling: tuple[int, int, int]
    cell: UnitCell
    origin: tuple[float, float, float]
    space_group: int
    labels: list[str]
    extended_header_bytes: int
    mode: int
    value_range: tuple[float, float]

@final
class MrcBlockReader:
    def __init__(self, path: str, *, options: MrcBlockOptions | None = ...) -> None: ...
    descriptor: MrcMapDescriptor
    def read_block(self, origin: tuple[int, int, int], dimensions: tuple[int, int, int]) -> NDArray[float32]: ...

@final
class MapBrickId:
    value: int
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __index__(self) -> int: ...

@final
class MrcBrickBudget:
    def __init__(self, *, max_payload_bytes: int = ..., max_working_set_bytes: int = ...) -> None: ...
    max_payload_bytes: int
    max_working_set_bytes: int

@final
class MrcBrickOptions:
    def __init__(self, *, interior: tuple[int, int, int] = ..., halo: int = ..., generation: int = ..., budget: MrcBrickBudget | None = ...) -> None: ...
    interior: tuple[int, int, int]
    halo: int
    generation: int
    budget: MrcBrickBudget

@final
class MapBrickAddress:
    origin: tuple[int, int, int]
    mip: int

@final
class MapBrickShape:
    stored: tuple[int, int, int]
    interior: tuple[int, int, int]
    halo: int
    voxel_count: int

@final
class ScalarBrickMetadata:
    id: MapBrickId
    address: MapBrickAddress
    shape: MapBrickShape
    minimum: float
    maximum: float
    generation: int

@final
class MrcBrickDescriptor:
    dataset: DatasetId
    chunk: ChunkId
    metadata: ScalarBrickMetadata

@final
class ScalarBrickPayload:
    descriptor: MrcBrickDescriptor
    actual_range: tuple[float, float]
    payload_bytes: int
    values: NDArray[float32]

class MrcBrickProvider:
    def __init__(self, path: str, dataset: DatasetId, first_chunk: ChunkId, first_brick: MapBrickId, *, options: MrcBrickOptions | None = ...) -> None: ...
    dataset: DatasetId
    logical_extent: tuple[int, int, int]
    brick_count: int
    def descriptor(self, id: MapBrickId) -> MrcBrickDescriptor: ...
    def read_brick(self, id: MapBrickId) -> ScalarBrickPayload: ...

DEFAULT_MRC_BLOCK_MEMORY_LIMIT_BYTES: int
DEFAULT_MRC_BRICK_PAYLOAD_BYTES: int
DEFAULT_MRC_BRICK_WORKING_SET_BYTES: int
def read_cube(data: bytes) -> CubeGrid: ...
def read_dx(data: bytes) -> DensityMap: ...

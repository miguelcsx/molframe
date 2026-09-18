from os import PathLike
from numpy import float64, uint32
from numpy.ma import MaskedArray
from numpy.typing import NDArray

class DmsParticle:
    def __init__(self, *, atomic_number: int | None = ..., component: int | None = ..., nonbonded_type: int | None = ..., mass: float | None = ..., charge: float | None = ..., residue_id: int | None = ..., residue_name: str | None = ..., chain: str | None = ..., segment: str | None = ..., name: str | None = ..., insertion: str | None = ..., formal_charge: float | None = ..., occupancy: float | None = ..., b_factor: float | None = ..., temperature_group: int | None = ..., energy_group: int | None = ..., ligand_group: int | None = ..., bias_group: int | None = ...) -> None: ...
    atomic_number: int | None
    component: int | None
    nonbonded_type: int | None
    mass: float | None
    charge: float | None
    residue_id: int | None
    residue_name: str | None
    chain: str | None
    segment: str | None
    name: str | None
    insertion: str | None
    formal_charge: float | None
    occupancy: float | None
    b_factor: float | None
    temperature_group: int | None
    energy_group: int | None
    ligand_group: int | None
    bias_group: int | None

class DmsVersion:
    def __init__(self, major: int, minor: int) -> None: ...
    major: int
    minor: int

class DmsBond:
    def __init__(self, p0: int, p1: int, order: float) -> None: ...
    p0: int
    p1: int
    order: float

class DmsCell:
    def __init__(self, vectors: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]) -> None: ...
    vectors: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]

class DmsFrame:
    def __init__(self, positions: list[tuple[float, float, float]], *, velocities: list[tuple[float, float, float] | None] | None = ..., cell: DmsCell | None = ...) -> None: ...
    positions: list[tuple[float, float, float]]
    velocities: list[tuple[float, float, float] | None]
    cell: DmsCell | None

class DmsTopology:
    def __init__(self, particles: list[DmsParticle], bonds: list[DmsBond]) -> None: ...
    particles: list[DmsParticle]
    bonds: list[DmsBond]

class DmsSystem:
    def __init__(self, *, version: DmsVersion | None = ..., topology: DmsTopology | None = ..., frame: DmsFrame | None = ...) -> None: ...
    @staticmethod
    def read(path: str | PathLike[str]) -> DmsSystem: ...
    def write(self, path: str | PathLike[str]) -> None: ...
    def __len__(self) -> int: ...
    version: tuple[int, int] | None
    version_record: DmsVersion | None
    particles: list[DmsParticle]
    topology: DmsTopology
    frame: DmsFrame
    positions: NDArray[float64]
    velocities: MaskedArray
    bond_indices: NDArray[uint32]
    bond_orders: NDArray[float64]
    cell: NDArray[float64] | None

def read_dms(path: str | PathLike[str]) -> DmsSystem: ...
def write_dms(path: str | PathLike[str], system: DmsSystem) -> None: ...

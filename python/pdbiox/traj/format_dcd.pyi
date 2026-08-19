from typing import final
from .._trajectory import Timestep

@final
class DcdEndian:
    Little: DcdEndian
    Big: DcdEndian

@final
class DcdHeader:
    def __init__(self, frame_count: int, atom_count: int, start_step: int, save_interval: int, delta_akma: float, fixed_atom_count: int, titles: list[str], endian: DcdEndian, charmm: bool, has_unit_cell: bool, has_fourth_dimension: bool) -> None: ...
    frame_count: int
    atom_count: int
    start_step: int
    save_interval: int
    delta_akma: float
    fixed_atom_count: int
    titles: list[str]
    endian: DcdEndian
    charmm: bool
    has_unit_cell: bool
    has_fourth_dimension: bool

@final
class DcdTrajectory:
    def __init__(self, header: DcdHeader, frames: list[Timestep]) -> None: ...
    header: DcdHeader
    frames: list[Timestep]

@final
class DcdWriteOptions:
    def __init__(self, *, endian: DcdEndian | None = ..., title: str | None = ..., start_step: int | None = ..., save_interval: int | None = ..., delta_akma: float | None = ...) -> None: ...
    endian: DcdEndian
    title: str
    start_step: int
    save_interval: int
    delta_akma: float

class DcdError(Exception): ...
def parse_dcd(bytes: bytes) -> DcdTrajectory: ...
def write_dcd(frames: list[Timestep], options: DcdWriteOptions) -> bytes: ...

__all__: list[str]

from __future__ import annotations
from os import PathLike
from typing import Literal, Protocol, final

class Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

@final
class Trajectory:
    @property
    def format(self) -> str: ...
    @property
    def positions(self) -> Array: ...
    @property
    def times(self) -> Array: ...
    @property
    def n_frames(self) -> int: ...
    @property
    def n_atoms(self) -> int: ...
    def __len__(self) -> int: ...

def read(
    path: str | PathLike[str],
    *,
    format: Literal["xtc", "trr", "dcd", "tng", "gro", "xyz", "lammps_dump", "netcdf"]
    | None = None,
) -> Trajectory: ...

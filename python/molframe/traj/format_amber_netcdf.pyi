from typing import final
from .._trajectory import Timestep

@final
class AmberNetcdfPrecision:
    Single: AmberNetcdfPrecision
    Double: AmberNetcdfPrecision

@final
class AmberNetcdfMetadata:
    precision: AmberNetcdfPrecision
    program: str | None
    program_version: str | None

@final
class AmberNetcdfTrajectory:
    frames: list[Timestep]
    metadata: AmberNetcdfMetadata

@final
class AmberNetcdfWriteOptions:
    def __init__(self, *, precision: AmberNetcdfPrecision | None = ..., program: str | None = ..., program_version: str | None = ...) -> None: ...
    precision: AmberNetcdfPrecision
    program: str
    program_version: str

class AmberNetcdfError(Exception): ...
def parse_amber_netcdf(bytes: bytes) -> list[Timestep]: ...
def parse_amber_netcdf_record(bytes: bytes) -> AmberNetcdfTrajectory: ...
def write_amber_netcdf(frames: list[Timestep], options: AmberNetcdfWriteOptions) -> bytes: ...

__all__: list[str]

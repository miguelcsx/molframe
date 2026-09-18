from typing import final
from .._trajectory import Timestep

@final
class TrrPrecision:
    Single: TrrPrecision
    Double: TrrPrecision

@final
class TrrTrajectory:
    frames: list[Timestep]
    steps: list[int]
    precision: list[TrrPrecision]

@final
class TrrWriteOptions:
    def __init__(self, precision: TrrPrecision | None = ...) -> None: ...
    precision: TrrPrecision

class TrrError(Exception): ...
def parse_trr(bytes: bytes) -> TrrTrajectory: ...
def write_trr(frames: list[Timestep], options: TrrWriteOptions) -> bytes: ...
def write_trr_with_precisions(frames: list[Timestep], precisions: list[TrrPrecision]) -> bytes: ...

__all__: list[str]

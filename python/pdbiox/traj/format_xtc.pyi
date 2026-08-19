from typing import final
from .._trajectory import Timestep

@final
class XtcTrajectory:
    frames: list[Timestep]
    steps: list[int]
    precision: list[float]

@final
class XtcWriteOptions:
    def __init__(self, precision: float | None = ...) -> None: ...
    precision: float

class XtcError(Exception): ...
def parse_xtc(bytes: bytes) -> XtcTrajectory: ...
def write_xtc(frames: list[Timestep], options: XtcWriteOptions) -> bytes: ...
def write_xtc_with_precisions(frames: list[Timestep], precisions: list[float]) -> bytes: ...

__all__: list[str]

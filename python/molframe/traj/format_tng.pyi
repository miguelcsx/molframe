from os import PathLike
from .._trajectory import Timestep

class TngCompression:
    @staticmethod
    def uncompressed() -> TngCompression: ...
    @staticmethod
    def lossless() -> TngCompression: ...
    @staticmethod
    def lossy(precision: float) -> TngCompression: ...
    kind: str
    precision: float | None
class TngError(Exception): ...
class TngWriteOptions:
    def __init__(self, *, distance_unit_exponent: int = ..., compression: TngCompression | None = ..., hashes: bool = ...) -> None: ...
    distance_unit_exponent: int
    compression: TngCompression
    hashes: bool
class TngTrajectory:
    frames: list[Timestep]
    steps: list[int]
    distance_unit_exponent: int
    compression_precision: float
    compression: TngCompression
def parse_tng(path: str | PathLike[str]) -> TngTrajectory: ...
def write_tng(path: str | PathLike[str], frames: list[Timestep], options: TngWriteOptions) -> None: ...

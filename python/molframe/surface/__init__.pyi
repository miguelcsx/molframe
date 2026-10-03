from typing import Protocol

from .. import ExecutionContext

class Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

def sasa(
    coordinates: Array,
    radii: Array,
    *,
    probe: float = 1.4,
    points: int = 960,
    context: ExecutionContext | None = None,
) -> Array: ...
def lee_richards(
    coordinates: Array,
    radii: Array,
    *,
    probe: float = 1.4,
    slices: int = 20,
    context: ExecutionContext | None = None,
) -> Array: ...
def cavities(
    coordinates: Array, radii: Array, *, probe: float = 1.4, resolution: float = 0.5
) -> list[tuple[float, tuple[float, float, float], int]]: ...

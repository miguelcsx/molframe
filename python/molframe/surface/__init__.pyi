from __future__ import annotations
from typing import Protocol

class Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

def sasa(
    coordinates: Array, radii: Array, *, probe: float = 1.4, points: int = 960
) -> Array: ...
def lee_richards(
    coordinates: Array, radii: Array, *, probe: float = 1.4, slices: int = 20
) -> Array: ...
def cavities(
    coordinates: Array, radii: Array, *, probe: float = 1.4, resolution: float = 0.5
) -> list[tuple[float, tuple[float, float, float], int]]: ...

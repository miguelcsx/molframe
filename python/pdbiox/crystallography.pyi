from typing import final
from numpy import float64
from numpy.typing import NDArray

@final
class UnitCell:
    def __init__(self, lengths: list[float], angles: list[float]) -> None: ...
    lengths: list[float]
    angles: list[float]
    fractional_to_cartesian_matrix: list[list[float]]
    cartesian_to_fractional_matrix: list[list[float]]
    def to_cartesian(self, coordinates: NDArray[float64]) -> NDArray[float64]: ...
    def to_fractional(self, coordinates: NDArray[float64]) -> NDArray[float64]: ...

@final
class SymmetryOperation:
    id: str
    rotation: list[list[int]]
    translation: list[tuple[int, int]]

@final
class SpaceGroup:
    hall_number: int
    international_number: int
    international_short: str
    international_full: str
    hall_symbol: str
    choice: str
    operations: list[SymmetryOperation]

def space_group_by_number(hall_number: int) -> SpaceGroup: ...
def space_group_by_symbol(hall_symbol: str) -> SpaceGroup: ...
def space_group_settings(international_number: int) -> list[SpaceGroup]: ...

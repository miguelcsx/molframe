from typing import Literal, Protocol

from .. import ExecutionContext, Float32Array, UInt32Array
from ..crystal import UnitCell

class Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

def cross_pairs(
    first: Array,
    second: Array,
    cutoff: float,
    *,
    cell: UnitCell | None = None,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    context: ExecutionContext | None = None,
) -> tuple[UInt32Array, UInt32Array, Float32Array]:
    """Return sorted cross-set pairs with indices local to each coordinate array.

    With `cell`, distances are minimum-image through that periodic cell.
    """

def neighbor_pairs(
    coordinates: Array,
    cutoff: float,
    *,
    cell: UnitCell | None = None,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    context: ExecutionContext | None = None,
) -> tuple[UInt32Array, UInt32Array, Float32Array]: ...

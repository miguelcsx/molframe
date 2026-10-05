from typing import Literal, Protocol

from .. import ExecutionContext, Float32Array, UInt32Array

class Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

def cross_pairs(
    first: Array,
    second: Array,
    cutoff: float,
    *,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    context: ExecutionContext | None = None,
) -> tuple[UInt32Array, UInt32Array, Float32Array]:
    """Return sorted cross-set pairs with indices local to each coordinate array."""

def neighbor_pairs(
    coordinates: Array,
    cutoff: float,
    *,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    context: ExecutionContext | None = None,
) -> tuple[UInt32Array, UInt32Array, Float32Array]: ...

from typing import Literal, Protocol

from .. import ExecutionContext

class Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

def neighbor_pairs(
    coordinates: Array,
    cutoff: float,
    *,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    context: ExecutionContext | None = None,
) -> tuple[Array, Array, Array]: ...

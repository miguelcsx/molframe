from typing import Literal, Protocol

class Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

def neighbor_pairs(
    coordinates: Array,
    cutoff: float,
    *,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
) -> tuple[Array, Array, Array]: ...

from collections.abc import Iterator
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

class PairChunks(Iterator[tuple[UInt32Array, UInt32Array, Float32Array]]):
    def __iter__(self) -> PairChunks: ...
    def __next__(self) -> tuple[UInt32Array, UInt32Array, Float32Array]: ...

def neighbor_pairs_chunks(
    coordinates: Array,
    cutoff: float,
    *,
    chunk_size: int = 1_048_576,
    cell: UnitCell | None = None,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    context: ExecutionContext | None = None,
) -> PairChunks:
    """Yield `(first, second, distance)` chunks of the sorted pair list.

    Concatenated, the chunks equal `neighbor_pairs`. Each chunk holds at most
    `chunk_size` pairs.
    """

def knn(
    points: Array,
    k: int,
    *,
    query: Array | None = None,
    max_distance: float | None = None,
    cell: UnitCell | None = None,
    exclude_self: bool | None = None,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    context: ExecutionContext | None = None,
) -> tuple[UInt32Array, Float32Array]:
    """Return `(indices, distances)` of the `k` nearest points, shape `(n_query, k)`.

    Rows ascend in distance, ties by index. Missing neighbours hold
    `4294967295` and `inf`. Without `query` the points query themselves and
    `exclude_self` defaults to true. Distances are float32 (minimum-image with
    `cell`); a target tied with the k-th within float32 rounding of the scipy
    float64 distance may be ranked either way. Without `max_distance` a periodic
    search covers two lattice shells of images.
    """

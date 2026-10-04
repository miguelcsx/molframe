from collections.abc import Mapping, Sequence
from os import PathLike
from typing import Literal, Protocol

from .. import ExecutionContext, Structure

class Float32Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class Int64Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class CoordinateTensor:
    """The coordinates of a structure, importable by ``from_dlpack``.

    ``numpy.from_dlpack``, ``torch.from_dlpack`` and ``jax.dlpack.from_dlpack`` accept
    this object. DLPack has no read-only flag, so every import copies the
    ``(atoms, 3)`` ``float32`` coordinates and the consumer owns its tensor.
    """

    @property
    def shape(self) -> tuple[int, int]: ...
    @property
    def dtype(self) -> str: ...
    def __dlpack__(
        self,
        *,
        stream: object | None = None,
        max_version: tuple[int, int] | None = None,
        dl_device: tuple[int, int] | None = None,
        copy: bool | None = None,
    ) -> object: ...
    def __dlpack_device__(self) -> tuple[int, int]: ...

class Graph:
    """Node and edge tensors in the layout PyG and DGL consume.

    ``edge_index`` is ``(2, edges)`` ``int64``, every source then every target.
    ``node_features`` and ``edge_features`` are ``float32`` with one named column
    per requested feature. Every array is the caller's own copy.
    """

    @property
    def node_count(self) -> int: ...
    @property
    def edge_count(self) -> int: ...
    @property
    def edge_index(self) -> Int64Array: ...
    @property
    def node_features(self) -> Float32Array: ...
    @property
    def edge_features(self) -> Float32Array: ...
    @property
    def node_feature_names(self) -> list[str]: ...
    @property
    def edge_feature_names(self) -> list[str]: ...

class ManifestEntry:
    """One structure of a dataset and what is known of it without loading it."""

    def __init__(
        self,
        id: str,
        path: str | PathLike[str],
        atom_count: int,
        *,
        resolution: float | None = None,
        method: str | None = None,
        deposition_date: str | None = None,
        sequence: str | None = None,
        structure_cluster: str | None = None,
        tags: Sequence[str] = ...,
        statistics: Mapping[str, float] = ...,
    ) -> None: ...
    @property
    def id(self) -> str: ...
    @property
    def path(self) -> str: ...
    @property
    def atom_count(self) -> int: ...
    @property
    def resolution(self) -> float | None: ...
    @property
    def method(self) -> str | None: ...
    @property
    def deposition_date(self) -> str | None: ...
    @property
    def sequence(self) -> str | None: ...
    @property
    def structure_cluster(self) -> str | None: ...
    @property
    def tags(self) -> list[str]: ...
    @property
    def statistics(self) -> dict[str, float]: ...

class Dataset:
    """An immutable, lazy list of structures described by a manifest.

    Filtering, splitting and batching read only the manifest; ``load`` reads one
    structure.
    """

    def __init__(self, entries: Sequence[ManifestEntry]) -> None: ...
    @staticmethod
    def from_manifest(path: str | PathLike[str]) -> Dataset: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> ManifestEntry: ...
    @property
    def entries(self) -> list[ManifestEntry]: ...
    def filter(
        self,
        *,
        resolution_below: float | None = None,
        method: str | None = None,
        minimum_atoms: int | None = None,
        maximum_atoms: int | None = None,
        tag: str | None = None,
    ) -> Dataset: ...
    def split(
        self,
        strategy: Literal["sequence_identity", "structural_cluster", "temporal", "random"],
        *,
        ratios: tuple[float, float, float],
        threshold: float | None = None,
        seed: int | None = None,
    ) -> DatasetSplit: ...
    def batches(self, size: int) -> list[Dataset]: ...
    def load(self, index: int) -> Structure: ...

class DatasetSplit:
    @property
    def train(self) -> Dataset: ...
    @property
    def validation(self) -> Dataset: ...
    @property
    def test(self) -> Dataset: ...
    @property
    def warnings(self) -> list[str]: ...

def coordinates(structure: Structure) -> CoordinateTensor: ...
def graph(
    structure: Structure,
    *,
    nodes: Literal["atoms", "residues"],
    edges: Literal["bonds", "contacts", "radius", "k_nearest"],
    direction: Literal["undirected", "symmetric", "directed"],
    cutoff: float | None = None,
    neighbors: int | None = None,
    node_features: Sequence[str] = ...,
    edge_features: Sequence[str] = ...,
    missing: float | None = None,
    backend: str = "auto",
    periodic: bool = False,
    context: ExecutionContext | None = None,
) -> Graph: ...
def write_atoms_ipc(
    structure: Structure,
    path: str | PathLike[str],
    *,
    metadata: Mapping[str, str] = ...,
) -> None: ...
def write_atoms_parquet(
    structure: Structure,
    path: str | PathLike[str],
    *,
    metadata: Mapping[str, str] = ...,
) -> None: ...

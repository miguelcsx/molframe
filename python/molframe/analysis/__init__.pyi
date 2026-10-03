from collections.abc import Sequence
from typing import Literal, Protocol

from .. import Analysis, AnalysisPolicy, ExecutionContext, Structure, Table

class ArrayColumn(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class ContactTable:
    def __len__(self) -> int: ...
    def __arrow_c_stream__(self, requested_schema: object | None = ...) -> object: ...
    @property
    def first(self) -> ArrayColumn: ...
    @property
    def second(self) -> ArrayColumn: ...
    @property
    def distance(self) -> ArrayColumn: ...

def atom_contacts(
    value: object,
    cutoff: float,
    *,
    backend: str = ...,
    context: ExecutionContext | None = None,
) -> ContactTable: ...
def dssp(structure: Structure) -> Table:
    """Native DSSP states as residue indices and stable SecondaryStructure integer codes."""

def contacts(
    structure: Structure,
    cutoff: float,
    *,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[ContactTable]: ...
def hydrogen_bonds(
    structure: Structure,
    *,
    max_distance: float = 3.5,
    min_angle: float = 120.0,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]: ...
def salt_bridges(
    structure: Structure,
    *,
    max_distance: float = 4.0,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]: ...

class GridSpec:
    def __init__(
        self, voxel_to_world: Sequence[Sequence[float]], dimensions: Sequence[int]
    ) -> None: ...
    @property
    def voxel_to_world(self) -> list[list[float]]: ...
    @property
    def dimensions(self) -> list[int]: ...

class ScalarGrid:
    @property
    def spec(self) -> GridSpec: ...
    @property
    def values(self) -> list[float]: ...

def contact_potential(
    structure: Structure, charges: Sequence[float], spec: GridSpec, *, cutoff: float = 12.0
) -> ScalarGrid: ...

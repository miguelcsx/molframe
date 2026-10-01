from __future__ import annotations
from typing import Literal, Protocol

from .. import Analysis, AnalysisPolicy, Structure, Table

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

def atom_contacts(value: object, cutoff: float, *, backend: str = ...) -> ContactTable: ...
def contacts(
    structure: Structure,
    cutoff: float,
    *,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
) -> Analysis[ContactTable]: ...
def hydrogen_bonds(
    structure: Structure,
    *,
    max_distance: float = 3.5,
    min_angle: float = 120.0,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
) -> Analysis[Table]: ...
def salt_bridges(
    structure: Structure,
    *,
    max_distance: float = 4.0,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
) -> Analysis[Table]: ...

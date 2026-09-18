from typing import final
from .. import Document
from .._io_types import MonomerLibraryReadError, RestraintError
from ..core.contract import Diagnostic
from .reflection import ReflectionTable

@final
class BondRestraint:
    atoms: tuple[str, str]
    target: float
    sigma: float
    kind: str | None

@final
class AngleRestraint:
    atoms: tuple[str, str, str]
    target: float
    sigma: float

@final
class TorsionRestraint:
    id: str | None
    atoms: tuple[str, str, str, str]
    target: float
    sigma: float
    period: int

@final
class PlaneAtomRestraint:
    atom: str
    sigma: float

@final
class PlaneRestraint:
    id: str
    atoms: list[PlaneAtomRestraint]

@final
class ChiralVolumeSign:
    Positive: ChiralVolumeSign
    Negative: ChiralVolumeSign
    Both: ChiralVolumeSign

@final
class ChiralRestraint:
    id: str | None
    atoms: tuple[str, str, str, str]
    sign: ChiralVolumeSign

@final
class MonomerRestraints:
    id: str
    atoms: list[str]
    bonds: list[BondRestraint]
    angles: list[AngleRestraint]
    torsions: list[TorsionRestraint]
    planes: list[PlaneRestraint]
    chirals: list[ChiralRestraint]

@final
class MonomerLibrary:
    def __init__(self) -> None: ...
    def get(self, id: str) -> MonomerRestraints | None: ...
    def __len__(self) -> int: ...
    is_empty: bool
    def items(self) -> list[tuple[str, MonomerRestraints]]: ...

def lower_monomer_library(document: Document) -> MonomerLibrary: ...
def read_monomer_library(data: bytes) -> tuple[MonomerLibrary, list[Diagnostic]]: ...
def lower_structure_factor_cif(document: Document) -> ReflectionTable: ...
def write_structure_factor_cif(table: ReflectionTable) -> str: ...

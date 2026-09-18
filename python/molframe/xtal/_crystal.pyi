from typing import final
from numpy import float64
from numpy.typing import NDArray
from . import Document

DEFAULT_CRYSTAL_IMAGE_LIMIT: int
DEFAULT_INSTANCE_LIMIT: int
ASSEMBLIES_EXTENSION: str
INSTANCE_ID_ANNOTATION: str
NCS_EXTENSION: str
SYMMETRY_EXTENSION: str

@final
class AffineTransform:
    def __init__(self, matrix: list[list[float]], translation: list[float]) -> None: ...
    @staticmethod
    def identity() -> AffineTransform: ...
    matrix: list[list[float]]
    translation: list[float]
    def is_finite(self) -> bool: ...
    def apply(self, position: tuple[float, float, float]) -> tuple[float, float, float]: ...
    def then(self, next: AffineTransform) -> AffineTransform: ...

@final
class Rational:
    def __init__(self, numerator: int, denominator: int) -> None: ...
    numerator: int
    denominator: int
    value: float

@final
class OperExpression:
    @staticmethod
    def parse(expression: str) -> OperExpression: ...
    @staticmethod
    def parse_with_limit(expression: str, limit: int) -> OperExpression: ...
    factors: list[list[str]]
    combination_count: int

@final
class Generator:
    oper_expression: OperExpression
    asym_ids: list[str]

@final
class Operator:
    id: str
    transform: object

@final
class AssemblyDef:
    id: str
    details: str | None
    method: str | None
    oligomeric: int | None
    generators: list[Generator]

@final
class AssemblySet:
    def __len__(self) -> int: ...
    is_empty: bool
    def get(self, id: str) -> AssemblyDef | None: ...
    def operator(self, id: str) -> Operator | None: ...
    def assemblies(self) -> list[AssemblyDef]: ...
    def operators(self) -> list[Operator]: ...

@final
class ChainInstance:
    source_chain: int
    transform: object
    instance_id: int

@final
class AtomInstance:
    source_atom: int
    transform: object
    instance_id: int
    def apply(self, position: tuple[float, float, float]) -> tuple[float, float, float]: ...

@final
class AssemblyNeighbor:
    first_instance: int
    first_atom: int
    second_instance: int
    second_atom: int
    distance_squared: float

@final
class AssemblyView:
    id: str
    instance_count: int
    def chains(self) -> list[ChainInstance]: ...
    def atoms(self) -> list[AtomInstance]: ...
    def positions(self, model: int = 0) -> list[tuple[AtomInstance, tuple[float, float, float] | None]]: ...
    def neighbors(self, model: int, cutoff: float, backend: object) -> list[AssemblyNeighbor]: ...

@final
class NcsCode:
    Given: NcsCode
    Generate: NcsCode

@final
class NcsOperator:
    id: str
    code: NcsCode
    details: str | None
    transform: AffineTransform

@final
class NcsSet:
    def __len__(self) -> int: ...
    is_empty: bool
    def get(self, id: str) -> NcsOperator | None: ...
    def operators(self) -> list[NcsOperator]: ...
    def generators(self) -> list[NcsOperator]: ...

@final
class NcsView:
    copy_count: int
    def atoms(self) -> list[tuple[int, NcsOperator]]: ...
    def positions(self, model: int = 0) -> list[tuple[int, NcsOperator, tuple[float, float, float] | None]]: ...

@final
class CrystalNeighbor:
    source_atom: int
    image_atom: int
    operation: int
    lattice: tuple[int, int, int]
    distance_squared: float

@final
class SymmetrySet:
    hall_number: int | None
    international_number: int | None
    hermann_mauguin: str | None
    hall: str | None
    crystal_system: str | None
    choice: str | None
    operations: list[SymmetryOperation]
    def is_empty(self) -> bool: ...

@final
class UnitCell:
    def __init__(self, lengths: list[float], angles: list[float]) -> None: ...
    lengths: list[float]
    angles: list[float]
    def is_placeholder(self) -> bool: ...
    fractional_to_cartesian_matrix: list[list[float]]
    cartesian_to_fractional_matrix: list[list[float]]
    def to_cartesian(self, coordinates: NDArray[float64]) -> NDArray[float64]: ...
    def to_fractional(self, coordinates: NDArray[float64]) -> NDArray[float64]: ...

@final
class CellTransform:
    def __init__(self, cell: UnitCell) -> None: ...
    fractional_to_cartesian_matrix: list[list[float]]
    cartesian_to_fractional_matrix: list[list[float]]
    def to_cartesian(self, fractional: tuple[float, float, float]) -> tuple[float, float, float]: ...
    def to_fractional(self, cartesian: tuple[float, float, float]) -> tuple[float, float, float]: ...

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

SpaceGroupSetting: type[SpaceGroup]

def space_group_by_number(hall_number: int) -> SpaceGroup: ...
def space_group_by_symbol(hall_symbol: str) -> SpaceGroup: ...
def space_group_by_hall(hall_symbol: str) -> SpaceGroup: ...
def space_group_setting(hall_number: int) -> SpaceGroup: ...
def space_group_settings(international_number: int) -> list[SpaceGroup]: ...
def lower_assemblies(document: Document) -> AssemblySet: ...
def lower_ncs(document: Document) -> NcsSet: ...
def lower_symmetry(document: Document) -> SymmetrySet: ...
def collect_crystal_neighbors(structure: object, symmetry: SymmetrySet, model: int = 0, cutoff: float = 5.0, *, limit: int | None = None, backend: object = ...) -> list[CrystalNeighbor]: ...

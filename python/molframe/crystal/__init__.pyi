from collections.abc import Mapping, Sequence
from os import PathLike
from typing import Literal, Protocol, Self, final

from .. import ExecutionContext, Structure

class IntArray(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class FloatArray(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class IntArray64(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class ComplexArray(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

@final
class UnitCell:
    def __init__(self, lengths: Sequence[float], angles: Sequence[float]) -> None: ...
    def reciprocal_vector(self, hkl: Sequence[int]) -> list[float]: ...
    def reciprocal_spacing_squared(self, hkl: Sequence[int]) -> float: ...
    def d_spacing(self, hkl: Sequence[int]) -> float: ...
    @property
    def lengths(self) -> tuple[float, float, float]: ...
    @property
    def angles(self) -> tuple[float, float, float]: ...
    def to_cartesian(self, fractional: Sequence[float]) -> list[float]: ...
    def to_fractional(self, cartesian: Sequence[float]) -> list[float]: ...

@final
class SymmetryOperation:
    """One operation on fractional coordinates: ``x' = W x + w``."""

    @property
    def expression(self) -> str:
        """The operation as written in the tables, such as ``1/2-x,1/2+y,-z``."""
    @property
    def rotation(self) -> tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]:
        """The integer rotation matrix ``W``, row by row."""
    @property
    def translation(self) -> tuple[float, float, float]:
        """The fractional translation ``w``."""
    @property
    def is_identity(self) -> bool: ...
    def apply(self, fractional: Sequence[float]) -> list[float]:
        """Return where the operation takes a fractional coordinate."""

@final
class SpaceGroup:
    def __init__(self, hall_number: int) -> None: ...
    @staticmethod
    def from_hermann_mauguin(symbol: str) -> SpaceGroup:
        """Return the group a Hermann-Mauguin symbol names, short or full.

        Case, spaces and screw-axis underscores do not matter, so a ``CRYST1`` spelling
        (``P 21 21 21``) resolves; a symbol of several settings gives the standard one.
        """
    @staticmethod
    def from_hall(symbol: str) -> SpaceGroup:
        """Return the setting a Hall symbol names."""
    @staticmethod
    def settings(international_number: int) -> list[SpaceGroup]:
        """Return every Hall setting of one International Tables type (1 to 230)."""
    def __len__(self) -> int: ...
    def reflection_symmetry(self, hkl: Sequence[int]) -> ReflectionSymmetry: ...
    @property
    def hall_symbol(self) -> str | None: ...
    @property
    def hall_number(self) -> int | None: ...
    @property
    def international_number(self) -> int | None: ...
    @property
    def hermann_mauguin(self) -> str | None: ...
    @property
    def choice(self) -> str | None:
        """The unique-axis, origin or cell choice that distinguishes settings of one type."""
    @property
    def operations(self) -> list[SymmetryOperation]: ...

@final
class ReflectionSymmetry:
    @property
    def centric(self) -> bool: ...
    @property
    def systematically_absent(self) -> bool: ...
    @property
    def epsilon_factor(self) -> int: ...

@final
class AssemblyInstance:
    @property
    def matrix(self) -> list[float]: ...
    @property
    def chains(self) -> list[str]: ...

@final
class AssemblyBond:
    @property
    def first_instance(self) -> int: ...
    @property
    def first_atom(self) -> int: ...
    @property
    def second_instance(self) -> int: ...
    @property
    def second_atom(self) -> int: ...
    @property
    def order(self) -> Literal["single"]: ...

def assembly_covalent_links(
    structure: Structure, id: str, *, model: int = 0, context: ExecutionContext | None = None
) -> list[AssemblyBond]: ...

@final
class ResolutionBins:
    def __new__(
        cls,
        cell: UnitCell,
        hkl: IntArray,
        *,
        bins: int = 20,
        method: Literal["equal_count", "dstar", "dstar2", "dstar3"] = "equal_count",
    ) -> Self: ...
    def __len__(self) -> int: ...
    @property
    def limits(self) -> FloatArray: ...
    @property
    def midpoints(self) -> FloatArray: ...
    def indices(self, inverse_d2: FloatArray) -> IntArray64: ...
    def d_min(self, bin: int) -> float: ...
    def d_max(self, bin: int) -> float: ...

@final
class ReducedCell:
    @property
    def lengths(self) -> tuple[float, float, float]: ...
    @property
    def angles(self) -> tuple[float, float, float]: ...
    @property
    def change_of_basis(self) -> list[list[int]]: ...
    @property
    def iterations(self) -> int: ...
    @property
    def converged(self) -> bool: ...

def reduce_cell(
    lengths: Sequence[float],
    angles: Sequence[float],
    *,
    epsilon: float = 1e-9,
    iteration_limit: int = 100,
) -> ReducedCell: ...
def normalizers(
    cell: UnitCell,
    space_group: SpaceGroup,
    hkl: IntArray,
    amplitudes: FloatArray,
    bins: ResolutionBins,
) -> FloatArray: ...
def assemblies(structure: Structure) -> list[str]: ...
def assembly(structure: Structure, id: str) -> list[AssemblyInstance]: ...
def structure_factors(structure: Structure, hkl: IntArray) -> ComplexArray: ...

__all__ = [
    "AssemblyInstance",
    "ReducedCell",
    "ReflectionSymmetry",
    "ResolutionBins",
    "SpaceGroup",
    "UnitCell",
    "assemblies",
    "assembly",
    "normalizers",
    "reduce_cell",
    "structure_factors",
]

class BoolArray3(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class UInt64Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class Float32Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

def space_group(structure: Structure) -> SpaceGroup | None:
    """Return the space group a structure carries, or ``None`` when it has none."""

@final
class MapStatistics:
    @property
    def count(self) -> int: ...
    @property
    def minimum(self) -> float: ...
    @property
    def maximum(self) -> float: ...
    @property
    def mean(self) -> float: ...
    @property
    def sigma(self) -> float:
        """The population standard deviation."""

@final
class DensityMap:
    """A scalar density grid: values on a regular grid of a unit cell."""

    def __init__(
        self,
        values: Float32Array,
        cell: UnitCell,
        *,
        sampling: Sequence[int] | None = None,
        starts: Sequence[int] = (0, 0, 0),
        origin: Sequence[float] = (0.0, 0.0, 0.0),
        space_group: int = 0,
        labels: Sequence[str] = (),
    ) -> None:
        """Build a map from ``values`` of shape ``(nz, ny, nx)`` (x varies fastest) on ``cell``.

        ``sampling`` is the number of grid intervals the cell is divided into along x, y and z
        (the grid's own size by default), ``starts`` the first stored grid point, ``origin`` the
        real-space origin in angstroms and ``space_group`` the International Tables number.
        """
    @property
    def dimensions(self) -> tuple[int, int, int]: ...
    @property
    def starts(self) -> tuple[int, int, int]: ...
    @property
    def sampling(self) -> tuple[int, int, int]: ...
    @property
    def cell(self) -> UnitCell: ...
    @property
    def origin(self) -> tuple[float, float, float]: ...
    @property
    def space_group(self) -> int: ...
    @property
    def labels(self) -> list[str]: ...
    @property
    def values(self) -> Float32Array:
        """The densities, ``(nz, ny, nx)``."""
    def statistics(self) -> MapStatistics:
        """Return the mean, population sigma, minimum and maximum over every voxel."""
    def masked_statistics(self, mask: BoolArray3) -> MapStatistics:
        """Return the same over the voxels where ``mask`` is true."""
    def histogram(self, bins: int, minimum: float, maximum: float) -> UInt64Array:
        """Return voxel counts in ``bins`` equal intervals; values outside the range are omitted."""
    def sample(
        self,
        positions: FloatArray,
        *,
        method: Literal["linear", "cubic"] = "linear",
        boundary: Literal["missing", "periodic"] = "missing",
    ) -> Float32Array:
        """Return the interpolated density at Cartesian positions; ``nan`` where there is none."""
    def to_mrc(self) -> bytes:
        """Return the map as MRC2014 bytes (mode 2, little-endian)."""

def read_mrc(source: str | bytes | PathLike[str]) -> DensityMap:
    """Read an MRC2014 or CCP4 scalar map; complex modes are refused, not truncated."""

@final
class ReflectionTable:
    """Reflections: Miller indices and the named columns measured at them."""

    def __init__(
        self,
        hkl: IntArray,
        columns: Mapping[str, tuple[str, Sequence[float]]],
        *,
        cell: UnitCell | None = None,
        space_group_number: int | None = None,
        space_group_name: str | None = None,
        title: str = "",
    ) -> None:
        """Build a table from Miller indices ``(n, 3)`` and ``(type, values)`` columns.

        The types are ``amplitude``, ``intensity``, ``standard_deviation``, ``phase``, ``flag``
        and ``real``; ``nan`` marks a value that was not recorded.
        """
    def __len__(self) -> int: ...
    @property
    def title(self) -> str: ...
    @property
    def cell(self) -> UnitCell | None: ...
    @property
    def space_group_number(self) -> int | None: ...
    @property
    def space_group_name(self) -> str | None: ...
    @property
    def row_count(self) -> int: ...
    @property
    def labels(self) -> list[str]: ...
    @property
    def types(self) -> dict[str, str]: ...
    @property
    def miller_indices(self) -> IntArray: ...
    @property
    def history(self) -> list[str]: ...
    def column(self, label: str) -> FloatArray:
        """Return one column as floats, ``nan`` where the value was not recorded."""
    def to_mtz(self) -> bytes: ...

def read_mtz(source: str | bytes | PathLike[str]) -> ReflectionTable:
    """Read an MTZ file."""

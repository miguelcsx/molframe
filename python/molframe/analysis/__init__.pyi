from collections.abc import Sequence
from typing import Literal, Protocol

from .. import Analysis, AnalysisPolicy, ExecutionContext, Float64Array, Query, Structure, Table
from .._structure import Selection

type _RadiusSet = Literal["bondi", "amber_united", "charmm", "alvarez"]

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
def contact_map(
    structure: Structure,
    *,
    cutoff: float,
    min_separation: int,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Residue pairs within ``cutoff`` Å, with their closest atom distance."""

def pi_stacking(
    structure: Structure,
    *,
    max_centre_distance: float,
    max_parallel_angle: float,
    min_perpendicular_angle: float,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Ring stacking; ``kind`` is 0 for parallel planes and 1 for T-shaped ones."""

def cation_pi(
    structure: Structure,
    *,
    max_distance: float,
    max_face_angle: float,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]: ...
def water_bridges(
    structure: Structure,
    *,
    max_distance: float,
    min_angle: float,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]: ...
def chain_interface(
    structure: Structure,
    *,
    first_chain: str,
    second_chain: str,
    cutoff: float,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Residues of two chains within ``cutoff`` Å (a name matches a label or author label)."""

def half_sphere_exposure(
    structure: Structure,
    *,
    radius: float,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]: ...
def nucleic_torsions(
    structure: Structure,
    *,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Seven torsions per nucleotide in degrees; ``NaN`` where one cannot be defined."""

def native_contacts(
    reference: Structure,
    target: Structure,
    *,
    cutoff: float,
    tolerance: float,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[dict[str, float]]: ...

class DensityGrid:
    """A weighted density on a rectilinear grid."""

    @property
    def origin(self) -> tuple[float, float, float]: ...
    @property
    def spacing(self) -> tuple[float, float, float]: ...
    @property
    def shape(self) -> tuple[int, int, int]: ...
    @property
    def excluded_weight(self) -> float:
        """Weight that fell outside the grid."""

    @property
    def values(self) -> Float64Array:
        """Weight per unit volume, shape ``(nx, ny, nz)``."""

def leaflets(
    structure: Structure,
    sites: str | Query | Selection,
    *,
    connection_distance: float,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Group the selected sites into connected components: ``site``, ``leaflet``."""

def radial_distribution(
    structure: Structure,
    first: str | Query | Selection,
    second: str | Query | Selection,
    *,
    minimum_distance: float,
    maximum_distance: float,
    bins: int,
    volume: float,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Pair counts and g(r) per shell: ``lower``, ``upper``, ``count``, ``distribution``."""

def coordination_numbers(
    structure: Structure,
    first: str | Query | Selection,
    second: str | Query | Selection,
    *,
    minimum_distance: float,
    maximum_distance: float,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Neighbours of each ``first`` atom among ``second`` in the shell: ``atom``, ``count``."""

def linear_density(
    structure: Structure,
    *,
    axis: Literal["x", "y", "z"],
    minimum: float,
    maximum: float,
    bins: int,
    weights: Literal["count", "mass"] | Sequence[float] | None = None,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Weighted profile along an axis: ``lower``, ``upper``, ``weight``, ``density``."""

def density_map(
    structure: Structure,
    *,
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float],
    shape: tuple[int, int, int],
    weights: Literal["count", "mass"] | Sequence[float] | None = None,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[DensityGrid]: ...
def pore_profile(
    structure: Structure,
    *,
    axis_origin: tuple[float, float, float],
    axis_direction: tuple[float, float, float],
    start: float,
    end: float,
    samples: int,
    search_radius: float,
    grid_spacing: float,
    probe_radius: float = 0.0,
    radii: _RadiusSet = "bondi",
    memory_limit_bytes: int = ...,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Largest probe fitting each slice: ``axial_coordinate``, ``centre_x/y/z``, ``radius``."""

def surface_contacts(
    structure: Structure,
    *,
    tolerance: float,
    probe: float,
    surface_density: float,
    minimum_area: float,
    radii: _RadiusSet = "bondi",
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Close atom pairs whose atoms are both exposed: ``first``, ``second``, ``distance``."""

def contacts_by_definition(
    structure: Structure,
    *,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    surface_tolerance: float = 0.2,
    surface_density: float = 4.0,
    surface_minimum_area: float = 0.25,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Return contacts under the policy's ``contact_def`` and ``vdw_radii``.

    Columns: ``first``, ``second``, ``distance``.

    ``distance:<tolerance>`` puts two atoms in contact when they are no further apart than
    their radii plus the tolerance; ``surface:<probe>`` when their expanded surfaces touch and
    both keep an exposed patch (the ``surface_*`` keywords say how that test samples).
    """

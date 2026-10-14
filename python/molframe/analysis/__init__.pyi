import os
from collections.abc import Sequence
from typing import Literal, Protocol, Self

from .. import (
    Analysis,
    AnalysisPolicy,
    ExecutionContext,
    Float32Array,
    Float64Array,
    Query,
    Structure,
    Table,
)
from .._structure import Selection

class HydrogenBondCandidateOptions:
    def __new__(
        cls,
        *,
        max_distance: float = 3.5,
        maximum_candidates: int = 2_000_000,
        backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    ) -> Self: ...

def hydrogen_bond_candidates(
    structure: Structure,
    start: Float32Array,
    end: Float32Array,
    *,
    options: HydrogenBondCandidateOptions | None = None,
    context: ExecutionContext | None = None,
) -> Table:
    """Conservative original-row candidates; requires instantaneous geometric filtering."""

def heavy_atom_hydrogen_bonds(
    structure: Structure,
    *,
    max_distance: float = 3.5,
    min_distance: float = 2.2,
    min_acceptor_angle: float | None = 90.0,
    covalent_distance: float = 1.9,
    include_same_residue: bool = False,
    include_sulfur: bool = False,
    role_table: str | None = None,
    role_file: str | os.PathLike[str] | None = None,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    context: ExecutionContext | None = None,
) -> Table:
    """Hydrogen bonds for structures without hydrogens, from heavy-atom geometry.

    Donors and acceptors are typed from residue and atom names (protein, water, nucleic
    acids; other groups use N/O or the CCD roles when annotated). A pair is kept when the
    donor-acceptor distance is within ``[min_distance, max_distance]``, the atoms are not
    1-2 or 1-3 covalent neighbours, and the donor-acceptor-antecedent angle is at least
    ``min_acceptor_angle`` (``None`` records it without filtering). Columns: ``donor``,
    ``acceptor``, ``distance`` and ``acceptor_angle`` (NaN without an antecedent), sorted
    by ``(donor, acceptor)``. No hydrogen is implied. Periodic images are not supported.

    ``role_table`` (text) or ``role_file`` (a path) supplies the typing table; with
    neither, the profile's ``heavy-roles`` entry selects one from the cache. Rust parses
    and validates it. ``Table.applied["roles"]`` names the table that applied (name,
    version, source, licence).
    """

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
    """Native DSSP states as residue indices and stable SecondaryStructure integer codes.

    Needs no role profile: amino-acid residues without explicit role columns take the
    built-in standard backbone roles (N, CA, C, O). Roles applied with
    ``chemistry.apply_polymer_role_profile`` are never overridden, and residues that match
    neither are left ``Unknown``.
    """

def contacts(
    structure: Structure,
    cutoff: float,
    *,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[ContactTable]:
    """Return atom pairs within ``cutoff``; indices are analysed atoms (``atom_origin``)."""

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
) -> Analysis[Table]:
    """Return oppositely charged group pairs within ``max_distance`` (N and O atoms).

    Charged groups are the titratable amino-acid groups (Asp/Glu carboxylate O, Lys NZ, Arg
    NE/NH1/NH2, chain termini) plus any atom with an annotated formal charge; histidine is
    not treated as charged. Bonded and same-residue pairs are excluded. An input with
    neither a titratable group nor an annotation is indeterminate, never an empty
    ``complete`` table. Indices refer to the analysed atoms (see ``atom_origin``).
    """

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
    """Residues of two chains within ``cutoff`` Å.

    A bare name matches a label or an author chain id and is refused when it means different
    chains in the two; write ``"label:A"`` or ``"auth:A"`` to choose the namespace.
    """

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
    tolerance: float = 0.0,
    match: Literal["index", "names"] = "index",
    mapping: Sequence[tuple[int, int]] | None = None,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[dict[str, float]]:
    """Return the fraction Q of the reference's contacts that ``target`` keeps.

    A reference pair within ``cutoff`` ångström is native; it is kept when the
    same pair is within ``cutoff + tolerance`` ångström in ``target``
    (``tolerance`` is additive, in ångström, not a multiplier). With
    ``tolerance=0.0`` a structure compared with itself has Q = 1. The contact is
    a hard cutoff over all atom pairs, not the smoothed Best-Hummer-Eaton
    switching function.

    ``match="index"`` pairs equal atom indices and requires equal atom counts.
    ``match="names"`` pairs atoms by chain, residue number, insertion code,
    residue name and atom name and counts only contacts between atoms present in
    both. ``mapping`` gives explicit ``(reference_atom, target_atom)`` index pairs
    and takes precedence over ``match``.

    The value holds ``native``, ``kept``, ``fraction`` and ``matched_atoms``.
    """

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
    radii: str | None = None,
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
    radii: str | None = None,
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

def sasa(
    structure: Structure,
    *,
    probe: float | None = None,
    method: str | None = None,
    n_points: int | None = None,
    hetero: str | None = None,
    metals: str | None = None,
    radii: str | None = None,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Float64Array]:
    """Return the solvent-accessible surface area of every analysed atom, in square angstroms.

    Every option left as ``None`` takes the Rust default shared with the Rust facade and the
    ``molframe sasa`` command: Shrake-Rupley, 960 points, a 1.4 angstrom probe, the polymer
    alone, and the van der Waals metal fallback.

    ``method`` is ``"shrake-rupley"`` or ``"lee-richards"``; ``n_points`` is the points per
    atom for the first and the slices per atom for the second (960 and 20 when ``None``).
    ``hetero`` chooses which heteroatoms are measured and occlude: ``"none"`` (the default,
    as NACCESS and FreeSASA), ``"solute"`` (ligands and ions), ``"all"``, or a comma-separated
    union of ``"waters"``, ``"ligands"`` and ``"ions"``. ``metals`` says where a metal the
    radius set lacks takes its radius: ``"none"``, ``"vdw"`` or ``"ionic"``.

    The radii are the policy's ``vdw_radii``, the hydrogens that count are its ``hydrogens``,
    and under ``crystal:<radius>`` the area is that of the atom among its symmetry mates. Use
    ``atom_origin`` to tie the result to the input.

    Atoms the options leave out, and atoms with no positive radius, are ``NaN`` and do not
    occlude their neighbours. An included atom whose element has no radius in the set (iron
    under ``bondi`` with ``metals="none"``) makes the result ``partial`` with a ``W6004``
    warning naming the element and ``vdw_radii``.

    ``radii`` names the radius table: a cache entry (``<name>.radii``) or a file path. ``None``
    takes the table the ``policy``'s ``vdw_radii`` names, found by that name in the table cache
    (absent: ``E7102``, never a quiet fall back to the profile); with no policy either it takes
    the profile's ``vdw-radii`` line. A ``radii`` that names another table than the policy
    raises ``PolicyError`` (``E7103``). Rust parses and validates it; a table that cannot be
    supplied raises ``E7102`` (absent) or ``E2201`` (malformed), which name how to supply it. The
    result records the table that applied, by name, version, source and licence, in its
    ``vdw_radii`` assumption and its ``radii_table`` parameter.
    """

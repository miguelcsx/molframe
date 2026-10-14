from os import PathLike
from typing import Literal, Protocol, TypedDict, final

from .. import Analysis, AnalysisPolicy, ExecutionContext, Selection, Structure, Table

class UInt32Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class Float32Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

@final
class ClashTable:
    def __len__(self) -> int: ...
    @property
    def first(self) -> UInt32Array: ...
    @property
    def second(self) -> UInt32Array: ...
    @property
    def overlap(self) -> Float32Array: ...

def clashes(
    structure: Structure,
    *,
    tolerance: float = 0.4,
    radii: str | None = None,
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    exclude_one_four: bool = True,
    perceive_bonds: bool = True,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[ClashTable]:
    """Pairs whose van der Waals spheres overlap by more than ``tolerance`` Å.

    1-2, 1-3 and (unless ``exclude_one_four`` is false) 1-4 bonded pairs, and
    pairs in different alternate locations, are never reported. ``perceive_bonds``
    completes the bond table with intra-residue bonds, peptide/phosphodiester
    links and disulfides before the exclusion, which a structure read from a
    coordinate file needs because it carries no bonds.
    """

def bond_length_deviations(
    structure: Structure,
    *,
    tolerance: float,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Bonds off their covalent expectation by more than ``tolerance`` Å."""

def cis_peptides(
    structure: Structure,
    *,
    threshold_degrees: float,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Peptide bonds with ``|omega|`` at most ``threshold_degrees``."""

def quality_flags(
    structure: Structure,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Atoms with a zero or out-of-range occupancy or a negative B factor (``issue`` 0, 1, 2)."""

def valence(
    structure: Structure,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Atoms with more bonds than their element allows."""

def planarity(
    structure: Structure,
    *,
    max_deviation: float,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[Table]:
    """Aromatic rings that deviate from their plane by more than ``max_deviation`` Å."""

def completeness(
    structure: Structure,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[list[dict[str, object]]]:
    """Per polymer chain, the modelled residues against the canonical sequence."""

def altloc_occupancy_sums(
    structure: Structure,
    *,
    expected_sum: float,
    tolerance: float,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[dict[str, object]]:
    """Alternate-location groups whose occupancies do not sum to ``expected_sum``."""

def ligand_geometry(
    structure: Structure,
    *,
    tolerance: float,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[dict[str, object]]:
    """Heterogen bonds off their covalent expectation by more than ``tolerance`` Å."""

def b_factor_distribution(
    structure: Structure,
    *,
    outlier_standard_deviations: float,
    selection: Selection | None = None,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[dict[str, object]]:
    """Summarise the recorded B factors and list the outliers."""

class Float32Matrix(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class TableInfo(TypedDict):
    name: str
    version: str
    source: str
    licence: str

@final
class ValidationReport:
    """Per-residue arrays, counts and applied-table provenance of one check.

    ``arrays`` values are numpy arrays that share memory with the report.
    ``tables`` lists the name, version, source and licence of every table or
    policy that decided the result.
    """

    @property
    def check(self) -> str: ...
    @property
    def arrays(self) -> dict[str, Float32Array | UInt32Array | Float32Matrix]: ...
    @property
    def summary(self) -> dict[str, float]: ...
    @property
    def tables(self) -> list[TableInfo]: ...
    def __getitem__(self, name: str) -> Float32Array | UInt32Array | Float32Matrix: ...

def protein_conformation(
    structure: Structure,
    *,
    rama_table: str | PathLike[str] | None = None,
    rotamer_table: str | PathLike[str] | None = None,
    chi_table: str | PathLike[str] | None = None,
    context: ExecutionContext | None = None,
) -> ValidationReport:
    """Ramachandran, rotamer and omega classes per residue.

    Classes are 0 outlier, 1 allowed, 2 favoured and 255 not evaluated; omega
    classes are 0 trans, 1 cis, 2 twisted. Each table is the
    ``*_table`` file you pass or the cache entry the profile names; none is built in,
    and a missing one raises with the code E7102 and how to supply it. Rust parses
    and validates the file. ``tables`` on the report lists what applied.
    """

def geometry_z(
    structure: Structure,
    *,
    cutoff: float = 4.0,
    geometry_table: str | PathLike[str] | None = None,
    chi_table: str | PathLike[str] | None = None,
    context: ExecutionContext | None = None,
) -> ValidationReport:
    """Bond and angle z-scores against the selected geometry table (file or profile entry)."""

def nucleic_conformation(
    structure: Structure,
    *,
    nucleic_table: str | PathLike[str] | None = None,
    context: ExecutionContext | None = None,
) -> ValidationReport:
    """Sugar pucker and backbone suite sanity per nucleotide.

    The sanity filter flags delta/pucker mismatch and gamma outside the
    staggered windows; it does not assign named suites.
    """

def clashscore(
    structure: Structure,
    *,
    tolerance: float = 0.4,
    radii: str | None = None,
    hydrogens: Literal["auto", "ignore"] = "auto",
    polar_margin: float = 0.6,
    context: ExecutionContext | None = None,
) -> ValidationReport:
    """Clashes per 1000 atoms; hydrogen-aware when the model has hydrogens."""

from typing import Literal, Protocol, final

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
    radii: Literal["bondi", "amber_united", "charmm", "alvarez"] = "bondi",
    backend: Literal["auto", "cell", "kd_tree", "brute_force"] = "auto",
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[ClashTable]: ...
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

from collections.abc import Sequence
from os import PathLike
from typing import Literal, Protocol

from .. import ExecutionContext, Structure
from ..sequence import Scoring

class Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

type _Names = Literal["label", "auth"]

class DockQ:
    """A docking score in ``(0, 1]`` and its components."""

    @property
    def fnat(self) -> float: ...
    @property
    def ligand_rmsd(self) -> float: ...
    @property
    def interface_rmsd(self) -> float: ...
    @property
    def score(self) -> float: ...

class MappedComparison:
    """What corresponded when chains and atoms were mapped before scoring."""

    @property
    def chains(self) -> list[tuple[str, str, float]]: ...
    @property
    def atoms(self) -> int: ...
    @property
    def swapped_residues(self) -> int: ...
    @property
    def status(self) -> Literal["complete", "partial", "ambiguous", "indeterminate", "unknown"]: ...
    @property
    def chain_assignments_tried(self) -> int: ...

def tm_score(model: Array, reference: Array) -> float: ...
def gdt_ts(model: Array, reference: Array) -> float: ...
def gdt_ha(model: Array, reference: Array) -> float: ...
def weighted_rmsd(model: Array, reference: Array, weights: Array) -> float: ...
def lddt(
    model: Array,
    reference: Array,
    *,
    inclusion_radius: float = 15.0,
    tolerances: Sequence[float] = (0.5, 1.0, 2.0, 4.0),
    minimum_reference_distance: float = 0.0,
    context: ExecutionContext | None = None,
) -> float:
    """Score how well a model preserves the reference's local distances, in ``[0, 1]``."""

def dockq(
    model: Structure,
    native: Structure,
    *,
    receptor: str,
    ligand: str,
    contact_distance: float = 5.0,
    ligand_scale: float = 8.5,
    interface_scale: float = 1.5,
    chain_names: _Names = "label",
) -> DockQ: ...
def qs_score(
    model: Structure,
    native: Structure,
    *,
    first_chain: str,
    second_chain: str,
    contact_distance: float = 5.0,
    chain_names: _Names = "label",
) -> float: ...
def mapped_dockq(
    model: Structure,
    native: Structure,
    *,
    receptor: str,
    ligand: str,
    components: str | PathLike[str],
    components_version: str,
    scoring: Scoring,
    min_identity: float,
    automorphism_limit: int,
    contact_distance: float = 5.0,
    ligand_scale: float = 8.5,
    interface_scale: float = 1.5,
    chain_names: _Names = "label",
) -> tuple[DockQ, MappedComparison]: ...
def mapped_qs_score(
    model: Structure,
    native: Structure,
    *,
    first_chain: str,
    second_chain: str,
    components: str | PathLike[str],
    components_version: str,
    scoring: Scoring,
    min_identity: float,
    automorphism_limit: int,
    contact_distance: float = 5.0,
    chain_names: _Names = "label",
) -> tuple[float, MappedComparison]: ...

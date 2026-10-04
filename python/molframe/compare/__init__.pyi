from collections.abc import Sequence
from os import PathLike
from typing import Literal, Protocol, final

from .. import ExecutionContext, Structure, Table
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

@final
class CeAlignment:
    """One combinatorial-extension correspondence between two sets of guide atoms."""

    def __len__(self) -> int: ...
    @property
    def reference_indices(self) -> Array:
        """Indices into the reference guide atoms."""
    @property
    def mobile_indices(self) -> Array:
        """The matching indices into the mobile guide atoms."""
    @property
    def fragment_count(self) -> int: ...
    @property
    def similarity(self) -> float:
        """The path similarity; values nearer zero are better."""
    @property
    def z_score(self) -> float | None:
        """The empirical significance estimate, when a calibration profile was chosen."""
    @property
    def rmsd(self) -> float:
        """RMSD of the correspondence after one shared rigid fit."""

def ce_align(
    reference: Array,
    mobile: Array,
    *,
    window_size: int = 8,
    max_gap: int = 30,
    max_paths: int = 20,
    fragment_threshold: float = -3.0,
    path_threshold: float = -4.0,
    significance: bool = True,
    memory_limit: int = ...,
) -> list[CeAlignment]:
    """Return combinatorial-extension alignments of two sets of guide atoms, best first.

    The defaults are the original CE search; ``significance=False`` reports no z-score.
    """

def contact_areas(
    positions: Array,
    radii: Array,
    residues: Array,
    *,
    probe: float = 1.4,
    density: float = 4.0,
    context: ExecutionContext | None = None,
) -> Table:
    """Return contact areas between residues. Columns: ``first``, ``second``, ``area``.

    ``residues`` gives each atom an identifier chosen so that structures to be compared share
    them.
    """

@final
class CadScore:
    """The contact-area difference score of a model against a reference."""

    @property
    def score(self) -> float: ...
    @property
    def reference_area(self) -> float: ...
    @property
    def lost_area(self) -> float: ...
    def contacts(self) -> Table:
        """Return every reference contact: first, second, reference, model and lost area."""
    def local(self) -> Table:
        """Return the score of each residue: residue, reference area, lost area, score."""

def cad_score(reference: Table, model: Table) -> CadScore:
    """Return the CAD score of a model's contact areas against a reference's."""

@final
class ContactSimilarity:
    @property
    def shared(self) -> int: ...
    @property
    def union(self) -> int: ...
    @property
    def jaccard(self) -> float:
        """``shared / union``; one when both maps are empty."""

def contact_similarity(first: Array, second: Array) -> ContactSimilarity:
    """Return the Jaccard overlap of two ``(n, 2)`` contact maps of residue pairs."""

@final
class ChainAssignment:
    @property
    def primary(self) -> list[tuple[str, str, float]]:
        """``(reference, target, identity)`` of the assignment of the greatest total identity."""
    @property
    def alternatives(self) -> list[tuple[str, str, float]]:
        """The pairs left out that still meet the identity threshold."""

def assign_chains(
    reference: Structure,
    target: Structure,
    *,
    components: str | PathLike[str],
    components_version: str,
    scoring: Scoring,
    min_identity: float,
    chain_names: _Names = "label",
) -> ChainAssignment:
    """Return the chain of ``target`` each chain of ``reference`` corresponds to, by sequence."""

def map_sequence_to_structure(
    query: str,
    structure: Structure,
    *,
    components: str | PathLike[str],
    components_version: str,
    scoring: Scoring,
    chain_names: _Names = "label",
) -> list[tuple[int, int]]:
    """Return ``(query_position, residue_index)`` for each query residue found in the structure."""

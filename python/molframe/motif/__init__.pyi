from os import PathLike
from typing import Literal, final

import numpy as np
from numpy.typing import NDArray

from .. import AnalysisPolicy, ExecutionContext, Structure

@final
class ConstraintResult:
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> float | None: ...
    @property
    def deviation(self) -> float | None: ...
    @property
    def satisfied(self) -> bool | None: ...

@final
class MotifEvaluation:
    @property
    def mapping_index(self) -> int: ...
    @property
    def profile(self) -> str: ...
    @property
    def verdict(self) -> Literal["pass", "fail", "indeterminate"]: ...
    @property
    def constraints(self) -> list[ConstraintResult]: ...

@final
class MotifReport:
    @property
    def mapping_ambiguous(self) -> bool: ...
    @property
    def evaluations(self) -> list[MotifEvaluation]: ...
    def __len__(self) -> int: ...

@final
class MotifMatches:
    """Matches as residue indices: a ``(matches, components)`` uint32 array."""

    @property
    def components(self) -> list[str]: ...
    @property
    def residues(self) -> NDArray[np.uint32]: ...
    def __len__(self) -> int: ...

@final
class MatchIterator:
    """Lazy chunks of at most ``chunk_size`` matches; memory stays at one chunk."""

    @property
    def components(self) -> list[str]: ...
    def __iter__(self) -> MatchIterator: ...
    def __next__(self) -> NDArray[np.uint32]: ...

def evaluate(
    structure: Structure,
    specification: str | PathLike[str],
    *,
    components: str | PathLike[str] | None = None,
    version: str = "unversioned",
    limits: tuple[int | None, int] = (None, 256),
    constrained: bool = False,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> MotifReport: ...
def count(
    structure: Structure,
    specification: str | PathLike[str],
    *,
    components: str | PathLike[str] | None = None,
    version: str = "unversioned",
    constrained: bool = True,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> int: ...
def matches(
    structure: Structure,
    specification: str | PathLike[str],
    *,
    components: str | PathLike[str] | None = None,
    version: str = "unversioned",
    constrained: bool = True,
    limit: int | None = None,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> MotifMatches: ...
def first(
    structure: Structure,
    specification: str | PathLike[str],
    n: int,
    *,
    components: str | PathLike[str] | None = None,
    version: str = "unversioned",
    constrained: bool = True,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> MotifMatches: ...
def iter_matches(
    structure: Structure,
    specification: str | PathLike[str],
    *,
    chunk_size: int = 65536,
    components: str | PathLike[str] | None = None,
    version: str = "unversioned",
    constrained: bool = True,
    policy: AnalysisPolicy | None = None,
    context: ExecutionContext | None = None,
) -> MatchIterator: ...

__all__ = [
    "ConstraintResult",
    "MatchIterator",
    "MotifEvaluation",
    "MotifMatches",
    "MotifReport",
    "count",
    "evaluate",
    "first",
    "iter_matches",
    "matches",
]

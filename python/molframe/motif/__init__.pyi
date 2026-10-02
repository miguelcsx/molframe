from os import PathLike
from typing import Literal, final

from .. import AnalysisPolicy, Structure

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

def evaluate(
    structure: Structure,
    specification: str | PathLike[str],
    *,
    components: str | PathLike[str] | None = None,
    version: str = "unversioned",
    limits: tuple[int, int] = (64, 256),
    policy: AnalysisPolicy | None = None,
) -> MotifReport: ...

__all__: list[str]

from collections.abc import Callable, Sequence
from typing import Literal, final

from .. import Analysis, AnalysisPolicy

type Metric = Literal[
    "set",
    "absolute",
    "relative",
    "rms",
    "correlation",
    "ranking",
    "graph",
    "graph-nodes",
    "graph-edges",
    "flip",
]

@final
class PolicySpace:
    """A baseline policy and the decisions to vary, before they become runs."""

    def __init__(self, baseline: AnalysisPolicy | None = None, *, max_runs: int = 4096) -> None: ...
    def vary(self, field: str, values: Sequence[str]) -> PolicySpace:
        """Also vary ``field`` over ``values``, written as in a policy (``"biological:1"``)."""

    @property
    def cost(self) -> int:
        """The exact number of runs the space expands to."""

    def plan(self) -> AuditPlan:
        """Expand into runs; contradictory combinations are refused."""

@final
class AuditPlan:
    @property
    def cost(self) -> int: ...
    @property
    def fields(self) -> list[str]: ...
    @property
    def policies(self) -> list[AnalysisPolicy]: ...
    @property
    def coordinates(self) -> list[list[int]]: ...
    def __len__(self) -> int: ...

@final
class Effect:
    """What one decision explains of the variation among runs."""

    @property
    def field(self) -> str: ...
    @property
    def mean_change(self) -> float:
        """Mean distance between runs that differ in this decision only."""

    @property
    def comparisons(self) -> int: ...
    @property
    def share(self) -> float:
        """Fraction of the total variation this decision explains alone."""

@final
class Interaction:
    """What two decisions explain together beyond their separate effects."""

    @property
    def first(self) -> str: ...
    @property
    def second(self) -> str: ...
    @property
    def share(self) -> float: ...

@final
class AuditResult:
    @property
    def runs(self) -> list[Analysis[object]]: ...
    @property
    def policies(self) -> list[AnalysisPolicy]: ...
    @property
    def metric(self) -> str: ...
    @property
    def read(self) -> list[str]:
        """The decisions the analysis recorded as applied."""

    @property
    def indeterminate(self) -> list[int]:
        """Runs with no defensible answer."""

    @property
    def indeterminate_fraction(self) -> float: ...
    @property
    def agreement_with_first(self) -> float | None: ...
    @property
    def effects(self) -> list[Effect] | None:
        """``None`` while any run has no answer."""

    @property
    def interactions(self) -> list[Interaction] | None: ...
    @property
    def higher_order(self) -> float | None: ...
    @property
    def total_variation(self) -> float | None: ...
    @property
    def mean_distance(self) -> float | None: ...
    @property
    def max_distance(self) -> float | None: ...

def run(
    plan: AuditPlan,
    analyse: Callable[[AnalysisPolicy], Analysis[object]],
    *,
    metric: Metric,
    project: Callable[[object], object] | None = None,
) -> AuditResult:
    """Run ``analyse`` under every policy of ``plan`` and measure how far the answers move.

    ``project`` reduces an analysis' value to what ``metric`` compares: items for ``set`` and
    ``ranking``, a number for ``absolute`` and ``relative``, numbers for ``rms`` and
    ``correlation``, a category for ``flip``, ``(nodes, edges)`` for the graph metrics. A plan
    that varies a decision the analysis never applied is refused after the first run.
    """

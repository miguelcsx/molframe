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
    def vary(
        self,
        field: str,
        values: Sequence[str],
        *,
        rationale: str = "",
        evidence: str = "",
    ) -> PolicySpace:
        """Also vary ``field`` over ``values``, written as in a policy (``"biological:1"``).

        ``rationale`` says why those alternatives are the defensible ones and ``evidence`` what
        supports that; both travel with the plan.
        """

    def forbid(self, when: tuple[str, str], then_not: tuple[str, str]) -> PolicySpace:
        """Forbid a universe holding both ``when`` and ``then_not`` (each ``(decision, value)``)."""

    @property
    def cost(self) -> int:
        """The exact number of runs the space expands to."""

    def plan(self, *, constrained: bool = False) -> AuditPlan:
        """Expand into runs.

        The plan is the whole product of the alternatives, so a combination that cannot run or
        is forbidden is refused. With ``constrained=True`` such combinations are dropped and
        counted instead; the plan is then not the whole product, the effect and interaction
        shares do not apply, and the Shapley shares carry the attribution.
        """

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
    @property
    def balanced(self) -> bool:
        """Whether the plan is the whole product of the alternatives."""

    @property
    def skipped(self) -> int:
        """How many combinations a constrained plan dropped."""

    @property
    def decisions(self) -> list[Decision]: ...
    def __len__(self) -> int: ...

@final
class Decision:
    """A varied decision: its class of uncertainty and the case for varying it."""

    @property
    def field(self) -> str: ...
    @property
    def uncertainty(self) -> Literal["structural", "interpretive", "algorithmic", "numerical"]: ...
    @property
    def rationale(self) -> str: ...
    @property
    def evidence(self) -> str: ...

@final
class Attribution:
    """A decision's, or a class's, Shapley share of the variation among runs."""

    @property
    def name(self) -> str: ...
    @property
    def share(self) -> float: ...

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
    def balanced(self) -> bool | None:
        """Whether every combination was run, so effect and interaction shares add up."""

    @property
    def shapley(self) -> list[Attribution] | None:
        """Each decision's Shapley share; they sum to one, balanced plan or not."""

    @property
    def by_class(self) -> list[Attribution] | None:
        """The same attribution over the classes of uncertainty."""

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
    project: Callable[[Analysis[object]], object] | None = None,
) -> AuditResult:
    """Run ``analyse`` under every policy of ``plan`` and measure how far the answers move.

    ``project`` reduces an ``Analysis`` (its value, and for instance its ``atom_origin``) to what
    ``metric`` compares: items for ``set`` and
    ``ranking``, a number for ``absolute`` and ``relative``, numbers for ``rms`` and
    ``correlation``, a category for ``flip``, ``(nodes, edges)`` for the graph metrics. A plan
    that varies a decision the analysis never applied is refused after the first run.
    """

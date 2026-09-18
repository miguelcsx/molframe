"""Typed declarative trajectory operations backed by the native facade."""

from typing import Any, final

from numpy import float32, float64, uint64
from numpy.typing import NDArray

from ._trajectory_ensemble import FrameAlignment
from .core.contract import Analysis
from .query import AnalysisPolicy


@final
class MeanSquaredDisplacementSeries:
    @property
    def lags(self) -> NDArray[Any]: ...

    @property
    def observations(self) -> NDArray[uint64]: ...

    @property
    def values(self) -> NDArray[float64]: ...


@final
class RmsdToReference:
    def __init__(
        self,
        *,
        frames: NDArray[float32],
        reference: int,
        alignment: FrameAlignment | None = None,
        policy: AnalysisPolicy | None = ...,
    ) -> None: ...

    def execute(self) -> Analysis[NDArray[float64]]: ...

    def explain(self) -> dict[str, Any]: ...

    def to_dict(self) -> dict[str, Any]: ...

    @staticmethod
    def from_dict(config: dict[str, Any]) -> RmsdToReference: ...

    def __repr__(self) -> str: ...


@final
class MeanSquaredDisplacementOp:
    def __init__(
        self,
        *,
        frames: NDArray[float32],
        maximum_lag: int,
        atoms: NDArray[Any] | None = ...,
        policy: AnalysisPolicy | None = ...,
    ) -> None: ...

    def execute(self) -> Analysis[MeanSquaredDisplacementSeries]: ...

    def explain(self) -> dict[str, Any]: ...

    def to_dict(self) -> dict[str, Any]: ...

    @staticmethod
    def from_dict(config: dict[str, Any]) -> MeanSquaredDisplacementOp: ...

    def __repr__(self) -> str: ...


@final
class PairwiseFittedRmsd:
    def __init__(
        self,
        *,
        frames: NDArray[float32],
        memory_limit: int | None = None,
        policy: AnalysisPolicy | None = ...,
    ) -> None: ...

    def execute(self) -> Analysis[NDArray[float64]]: ...

    def explain(self) -> dict[str, Any]: ...

    def to_dict(self) -> dict[str, Any]: ...

    @staticmethod
    def from_dict(config: dict[str, Any]) -> PairwiseFittedRmsd: ...

    def __repr__(self) -> str: ...


@final
class GeneralizedProcrustesMean:
    def __init__(
        self,
        *,
        frames: NDArray[float32],
        tolerance: float | None = None,
        maximum_iterations: int | None = None,
        policy: AnalysisPolicy | None = ...,
    ) -> None: ...

    def execute(self) -> Analysis[NDArray[float32]]: ...

    def explain(self) -> dict[str, Any]: ...

    def to_dict(self) -> dict[str, Any]: ...

    @staticmethod
    def from_dict(config: dict[str, Any]) -> GeneralizedProcrustesMean: ...

    def __repr__(self) -> str: ...

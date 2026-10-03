from typing import Literal, Protocol, final

from .. import Analysis, AnalysisPolicy, ExecutionContext, Structure

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

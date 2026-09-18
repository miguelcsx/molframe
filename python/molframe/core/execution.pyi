from os import PathLike
from typing import final

DEFAULT_MEMORY_BUDGET_BYTES: int

@final
class MemoryBudget:
    def __init__(self, bytes: int = ...) -> None: ...
    bytes: int

@final
class CancellationToken:
    def __init__(self) -> None: ...
    is_cancelled: bool
    def cancel(self) -> None: ...
    def check(self) -> None: ...

@final
class ScratchPolicy:
    def __init__(self, max_bytes: int = ...) -> None: ...
    max_bytes: int

@final
class TempStoragePolicy:
    def __init__(self) -> None: ...
    @staticmethod
    def disabled() -> TempStoragePolicy: ...
    @staticmethod
    def directory(root: str | PathLike[str], max_bytes: int) -> TempStoragePolicy: ...
    root: str | None
    max_bytes: int

@final
class MemoryLease:
    bytes: int

@final
class SpillFile:
    def append(self, payload: bytes) -> None: ...
    def finish(self) -> SpillArtifact: ...

@final
class SpillArtifact:
    bytes: int
    records: int
    path: str
    def reader(self) -> SpillReader: ...

@final
class SpillReader:
    def read_next(self, max_bytes: int = ...) -> bytes | None: ...

@final
class ExecutionContext:
    def __init__(
        self,
        *,
        memory_budget: MemoryBudget | None = ...,
        cancellation: CancellationToken | None = ...,
        scratch_policy: ScratchPolicy | None = ...,
        temp_storage_policy: TempStoragePolicy | None = ...,
        worker_budget: int | None = ...,
    ) -> None: ...
    memory_budget: MemoryBudget
    worker_budget: int
    cancellation: CancellationToken
    scratch_policy: ScratchPolicy
    temp_storage_policy: TempStoragePolicy
    spill_bytes: int
    total_spilled_bytes: int
    reserved_bytes: int
    peak_reserved_bytes: int
    live_batches: int
    peak_live_batches: int
    def reserve(self, bytes: int) -> MemoryLease: ...
    def create_spill_file(self, key: str) -> SpillFile: ...

@final
class WindowedFile:
    def __init__(self, path: str, max_window_bytes: int, context: ExecutionContext) -> None: ...
    capacity: int
    length: int
    def window(self, start: int, length: int) -> bytes: ...

@final
class BatchDemand:
    def __init__(self, max_rows: int, max_bytes: int) -> None: ...
    max_rows: int
    max_bytes: int
    can_accept_work: bool
    def accepts(self, rows: int, retained_bytes: int) -> bool: ...

@final
class Backpressure:
    Ready: Backpressure
    Pending: Backpressure
    Finished: Backpressure

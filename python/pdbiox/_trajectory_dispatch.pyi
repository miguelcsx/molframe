from ._native import UnitCell
from os import PathLike
from ._trajectory import Trajectory, TrajectoryFormat, TrajectoryUnits, TrajectoryWriteOptions, AmberRestartLayout
from ._trajectory_format_models import Timestep
from ._trajectory_runtime import RandomAccess, Units

class AmberAsciiReadOptions:
    def __init__(self, atom_count: int, *, periodic_box: bool = False) -> None: ...
    atom_count: int
    periodic_box: bool

class TrajectoryReadOptions:
    def __init__(
        self,
        *,
        format: TrajectoryFormat | None = None,
        length_to_angstrom: float | None = None,
        particle_group: str | None = None,
        units: TrajectoryUnits | None = None,
        amber_restart_layout: AmberRestartLayout | None = None,
        amber_ascii: AmberAsciiReadOptions | None = None,
    ) -> None: ...
    @staticmethod
    def defaults() -> TrajectoryReadOptions: ...

class TrajectoryReaderOptions:
    def __init__(
        self,
        *,
        format: TrajectoryFormat | None = None,
        memory_limit_bytes: int = ...,
    ) -> None: ...
    format: TrajectoryFormat | None
    memory_limit_bytes: int

class TrajectoryStreamReader:
    def __init__(self, path: str | PathLike[str], *, memory_limit_bytes: int = ...) -> None: ...
    format: str
    n_atoms: int
    n_frames: int | None
    units: Units
    random_access: RandomAccess
    def read_next(self) -> Timestep | None: ...
    def seek(self, frame: int) -> None: ...
    def rewind(self) -> None: ...

def read_trajectory(
    path: str | PathLike[str],
    *,
    options: TrajectoryReaderOptions | None = None,
) -> TrajectoryStreamReader: ...

def read_trajectory_materialized(
    path: str | PathLike[str],
    *,
    options: TrajectoryReadOptions | None = None,
) -> Trajectory: ...

def write_trajectory(
    path: str | PathLike[str],
    trajectory: Trajectory,
    *,
    options: TrajectoryWriteOptions | None = None,
) -> None: ...

from collections.abc import Callable
import numpy as np
import numpy.typing as npt
from .core.execution import ExecutionContext

def rmsf_stream(path: str | PathLike[str], context: ExecutionContext, *, frame_workspace_bytes: int = 8_000_000) -> npt.NDArray[np.float64]: ...
def rmsd_stream(path: str | PathLike[str], reference: npt.NDArray[np.float32], emit: Callable[[int, float | None, float], None], context: ExecutionContext, *, fit: bool = True, frame_workspace_bytes: int = 8_000_000) -> int: ...

class StreamFrame:
    frame: int
    time: float | None
    dt: float | None
    positions: npt.NDArray[np.float32]
    velocities: npt.NDArray[np.float32] | None
    forces: npt.NDArray[np.float32] | None
    cell: UnitCell | None

def run_analysis_stream(path: str | PathLike[str], consume: Callable[[StreamFrame], None], context: ExecutionContext, *, frame_workspace_bytes: int = 8_000_000) -> int: ...

from .spatial import SpatialSearchOptions

def contact_counts_stream(path: str | PathLike[str], cutoff: float, emit: Callable[[int, float | None, int], None], context: ExecutionContext, *, options: SpatialSearchOptions | None = None, frame_workspace_bytes: int = 8_000_000) -> int: ...

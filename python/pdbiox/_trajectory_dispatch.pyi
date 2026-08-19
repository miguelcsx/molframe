from os import PathLike
from ._trajectory import Trajectory, TrajectoryFormat, TrajectoryUnits, TrajectoryWriteOptions, AmberRestartLayout

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

def read_trajectory(
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

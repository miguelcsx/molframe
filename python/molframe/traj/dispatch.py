"""Path-based trajectory dispatch backed by the native Rust dispatcher."""

from .._native.traj import (
    AmberAsciiReadOptions,
    TrajectoryReaderOptions,
    TrajectoryStreamReader,
    TrajectoryReadOptions,
    TrajectoryWriteOptions,
    read_trajectory,
    read_trajectory_materialized,
    write_trajectory,
)

__all__ = [
    "AmberAsciiReadOptions",
    "TrajectoryReaderOptions",
    "TrajectoryStreamReader",
    "TrajectoryReadOptions",
    "TrajectoryWriteOptions",
    "read_trajectory",
    "read_trajectory_materialized",
    "write_trajectory",
]

"""Trajectory format dispatch backed by native Rust I/O."""

from .._trajectory import (
    AmberAsciiReadOptions,
    FormatMetadata,
    TrajectoryData,
    TrajectoryFormat,
    TrajectoryIoError,
    TrajectoryMetadata,
    TrajectoryReadOptions,
    TrajectoryWriteOptions,
    TrzWriteOptions,
    read_trajectory_materialized,
    write_trajectory,
)

__all__ = [
    "AmberAsciiReadOptions",
    "FormatMetadata",
    "TrajectoryData",
    "TrajectoryFormat",
    "TrajectoryIoError",
    "TrajectoryMetadata",
    "TrajectoryReadOptions",
    "TrajectoryWriteOptions",
    "TrzWriteOptions",
    "read_trajectory_materialized",
    "write_trajectory",
]

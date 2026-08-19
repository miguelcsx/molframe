"""Path-based trajectory dispatch backed by the native Rust dispatcher."""

from .._trajectory import (
    AmberAsciiReadOptions,
    TrajectoryReadOptions,
    TrajectoryWriteOptions,
    read_trajectory,
    write_trajectory,
)

__all__ = [
    "AmberAsciiReadOptions",
    "TrajectoryReadOptions",
    "TrajectoryWriteOptions",
    "read_trajectory",
    "write_trajectory",
]

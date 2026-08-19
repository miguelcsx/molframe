"""GROMACS TRR trajectory bindings."""

from .._trajectory import (
    TrrError,
    TrrPrecision,
    TrrTrajectory,
    TrrWriteOptions,
    parse_trr,
    write_trr,
    write_trr_with_precisions,
)

__all__ = [
    "TrrError",
    "TrrPrecision",
    "TrrTrajectory",
    "TrrWriteOptions",
    "parse_trr",
    "write_trr",
    "write_trr_with_precisions",
]

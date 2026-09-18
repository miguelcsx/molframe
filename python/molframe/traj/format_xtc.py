"""GROMACS XTC trajectory bindings."""

from .._trajectory import (
    XtcError,
    XtcTrajectory,
    XtcWriteOptions,
    parse_xtc,
    write_xtc,
    write_xtc_with_precisions,
)

__all__ = [
    "XtcError",
    "XtcTrajectory",
    "XtcWriteOptions",
    "parse_xtc",
    "write_xtc",
    "write_xtc_with_precisions",
]

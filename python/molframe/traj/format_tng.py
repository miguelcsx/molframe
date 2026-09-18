"""Trajectory Next Generation container format."""

from .._trajectory import TngCompression, TngError, TngTrajectory, TngWriteOptions, parse_tng, write_tng

__all__ = ["TngCompression", "TngError", "TngTrajectory", "TngWriteOptions", "parse_tng", "write_tng"]

"""Grace/xmgrace time-series parsing backed by the native Rust parser."""

from .._native.traj import XvgSeries, read_xvg

__all__ = ["XvgSeries", "read_xvg"]

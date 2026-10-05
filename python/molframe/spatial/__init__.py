"""Spatial search operations."""

from .._native import spatial as _native

neighbor_pairs = _native.neighbor_pairs
cross_pairs = _native.cross_pairs

__all__ = ["cross_pairs", "neighbor_pairs"]

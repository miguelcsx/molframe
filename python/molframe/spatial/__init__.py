"""Spatial search operations."""

from .._native import spatial as _native

neighbor_pairs = _native.neighbor_pairs

__all__ = ["neighbor_pairs"]

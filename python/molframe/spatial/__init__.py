"""Spatial search operations."""

from .._native import spatial as _native

neighbor_pairs = _native.neighbor_pairs
neighbor_pairs_chunks = _native.neighbor_pairs_chunks
cross_pairs = _native.cross_pairs
knn = _native.knn
PairChunks = _native.PairChunks

__all__ = ["PairChunks", "cross_pairs", "knn", "neighbor_pairs", "neighbor_pairs_chunks"]

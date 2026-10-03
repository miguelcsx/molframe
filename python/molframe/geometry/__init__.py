"""Numeric geometry kernels."""

from .._native import geometry as _native

centroid = _native.centroid
distance_matrix = _native.distance_matrix
rmsd = _native.rmsd

__all__ = ["centroid", "distance_matrix", "rmsd"]

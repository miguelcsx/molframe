"""Ensemble geometry backed by zero-copy native Rust array kernels."""

from .._trajectory import (
    DEFAULT_PAIRWISE_MEMORY_LIMIT,
    DEFAULT_PROCRUSTES_MAXIMUM_ITERATIONS,
    DEFAULT_PROCRUSTES_TOLERANCE,
    EnsembleDistanceMatrix,
    EnsembleGeometryError,
    FrameAlignment,
    generalized_procrustes_mean,
    pairwise_fitted_rmsd,
    pairwise_torus_distance,
    rmsd_to_reference,
)

__all__ = [
    "DEFAULT_PAIRWISE_MEMORY_LIMIT",
    "DEFAULT_PROCRUSTES_MAXIMUM_ITERATIONS",
    "DEFAULT_PROCRUSTES_TOLERANCE",
    "EnsembleDistanceMatrix",
    "EnsembleGeometryError",
    "FrameAlignment",
    "generalized_procrustes_mean",
    "pairwise_fitted_rmsd",
    "pairwise_torus_distance",
    "rmsd_to_reference",
]

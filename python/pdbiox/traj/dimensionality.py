"""Trajectory dimensionality reduction backed by native Rust kernels."""

from .._trajectory import (
    CartesianFit,
    PcaResult,
    analyse_cartesian_pca,
    analyse_dihedral_pca,
    cartesian_pca,
    torsion_pca,
)

__all__ = [
    "CartesianFit",
    "PcaResult",
    "analyse_cartesian_pca",
    "analyse_dihedral_pca",
    "cartesian_pca",
    "torsion_pca",
]

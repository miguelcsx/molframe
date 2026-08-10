"""Biopython-shaped entry points backed by native pdbiox code."""

from .._native import PDBParser

__all__ = ["PDBParser"]

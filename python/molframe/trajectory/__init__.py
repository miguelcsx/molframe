"""Trajectory operations."""
from .._native import trajectory as _native

Trajectory = _native.Trajectory
read = _native.read

__all__ = ["Trajectory", "read"]

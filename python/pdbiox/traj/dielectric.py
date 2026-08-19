"""Trajectory dielectric kernels backed by native Rust execution."""

from .._trajectory import DielectricEstimate, TrajectoryDielectricOptions as DielectricOptions, dielectric_from_dipoles

__all__ = ["DielectricEstimate", "DielectricOptions", "dielectric_from_dipoles"]

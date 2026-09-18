"""Diffusion-map embeddings backed by native Rust execution."""

from .._trajectory import DiffusionMap, analyse_diffusion_map, diffusion_map

__all__ = ["DiffusionMap", "analyse_diffusion_map", "diffusion_map"]

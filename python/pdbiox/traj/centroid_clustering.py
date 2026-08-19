"""Centroid clustering backed by the native Rust k-means kernel."""

from .._trajectory import KMeans, KMeansError, KMeansOptions, kmeans

__all__ = ["KMeans", "KMeansError", "KMeansOptions", "kmeans"]

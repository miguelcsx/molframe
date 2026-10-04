"""Interoperability: Arrow streams, DLPack coordinates, graph tensors and datasets."""

from .._native import interop as _native

CoordinateTensor = _native.CoordinateTensor
Dataset = _native.Dataset
DatasetSplit = _native.DatasetSplit
Graph = _native.Graph
ManifestEntry = _native.ManifestEntry
coordinates = _native.coordinates
graph = _native.graph
write_atoms_ipc = _native.write_atoms_ipc
write_atoms_parquet = _native.write_atoms_parquet

__all__ = [
    "CoordinateTensor",
    "Dataset",
    "DatasetSplit",
    "Graph",
    "ManifestEntry",
    "coordinates",
    "graph",
    "write_atoms_ipc",
    "write_atoms_parquet",
]

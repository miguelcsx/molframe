"""Lazy datasets and leakage-aware splits."""

from .._native import (
    ArrowStream, AtomArrowTable, AtomArrowTable as AtomTable,
    BondArrowTable, BondArrowTable as BondTable,
    ChainArrowTable, ChainArrowTable as ChainTable,
    ResidueArrowTable, ResidueArrowTable as ResidueTable,
    DLDataType, DLDevice, DLManagedTensor, DLTensor, Dataset, DatasetEntry,
    DatasetError, DatasetFilter,
    DatasetSplit, DatasetWarning, DlpackError, DlpackTensor, ExportCost, GraphError,
    LoadError, ManifestEntry, MolframeExtension, SplitOptions, SplitRatios,
    SplitStrategy, TableFileError, extension_name, graph as build_graph, write_atom_ipc,
    write_atom_ipc_with_metadata, write_atom_parquet, write_atom_parquet_with_metadata,
)
from .graph import (
    EdgeDirection, EdgeFeature, EdgeKind, Graph, GraphOptions, MissingFeaturePolicy,
    NodeFeature, NodeLevel, SpatialBackend,
)
from . import graph

__all__ = [
    "ArrowStream", "AtomTable", "AtomArrowTable", "BondTable", "BondArrowTable", "ChainTable", "ChainArrowTable", "ResidueTable", "ResidueArrowTable", "DLDataType", "DLDevice",
    "DLManagedTensor", "DLTensor", "Dataset", "DatasetEntry", "DatasetError",
    "DatasetFilter", "DatasetSplit", "DatasetWarning", "DlpackError", "DlpackTensor",
    "ExportCost", "GraphError", "LoadError", "ManifestEntry", "MolframeExtension",
    "EdgeDirection", "EdgeFeature", "EdgeKind", "Graph", "GraphOptions",
    "MissingFeaturePolicy", "NodeFeature", "NodeLevel", "SpatialBackend",
    "SplitOptions", "SplitRatios", "SplitStrategy", "TableFileError",
    "extension_name", "graph", "build_graph", "write_atom_ipc", "write_atom_ipc_with_metadata",
    "write_atom_parquet", "write_atom_parquet_with_metadata",
]

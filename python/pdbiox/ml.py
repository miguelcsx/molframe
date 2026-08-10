"""Lazy datasets and leakage-aware splits."""

from ._native import (
    Dataset, DatasetEntry, DatasetSplit, DatasetWarning, SplitStrategy,
    write_atom_ipc, write_atom_parquet,
)

__all__ = [
    "Dataset", "DatasetEntry", "DatasetSplit", "DatasetWarning", "SplitStrategy",
    "write_atom_ipc", "write_atom_parquet",
]

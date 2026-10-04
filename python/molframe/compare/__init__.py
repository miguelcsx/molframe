"""Structure comparison operations."""

from .._native import compare as _native

DockQ = _native.DockQ
MappedComparison = _native.MappedComparison
dockq = _native.dockq
gdt_ha = _native.gdt_ha
gdt_ts = _native.gdt_ts
lddt = _native.lddt
mapped_dockq = _native.mapped_dockq
mapped_qs_score = _native.mapped_qs_score
qs_score = _native.qs_score
tm_score = _native.tm_score
weighted_rmsd = _native.weighted_rmsd

__all__ = [
    "DockQ",
    "MappedComparison",
    "dockq",
    "gdt_ha",
    "gdt_ts",
    "lddt",
    "mapped_dockq",
    "mapped_qs_score",
    "qs_score",
    "tm_score",
    "weighted_rmsd",
]

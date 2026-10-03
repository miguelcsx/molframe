"""Structure comparison operations."""

from .._native import compare as _native

gdt_ha = _native.gdt_ha
gdt_ts = _native.gdt_ts
tm_score = _native.tm_score
weighted_rmsd = _native.weighted_rmsd

__all__ = ["gdt_ha", "gdt_ts", "tm_score", "weighted_rmsd"]

"""Structure validation operations."""

from .._native import validation as _native

ClashTable = _native.ClashTable
clashes = _native.clashes

__all__ = ["ClashTable", "clashes"]

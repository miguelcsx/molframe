"""Query completion at a source cursor."""

from .._native import query as _native

complete = _native.complete

__all__ = ["complete"]

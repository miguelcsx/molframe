"""UPGMA hierarchical clustering.

Rust names this module ``upgma``. Python retains ``seq.upgma`` as the canonical
callable, so the module surface lives at this unambiguous path.
"""

from . import upgma

__all__ = ["upgma"]

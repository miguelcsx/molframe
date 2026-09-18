"""Lee-Richards slice integration.

Rust names this module ``lee_richards``. Python retains
``surface.lee_richards`` as the canonical callable, so its module surface lives
at this unambiguous path.
"""

from . import lee_richards

__all__ = ["lee_richards"]

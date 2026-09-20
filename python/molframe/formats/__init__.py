"""Structure format namespaces."""
from .._native import formats as _native
from . import bcif, cif, modelcif, pdb
__all__ = ["bcif", "cif", "modelcif", "pdb"]

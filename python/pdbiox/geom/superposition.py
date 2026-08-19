"""Rigid least-squares superposition.

Rust names this module ``superpose``. Python retains ``geom.superpose`` as the
canonical callable, so the module surface lives at this unambiguous path.
"""

from . import SuperposeError, SuperposeOptions, Superposition, rmsd, rmsd_flat, superpose, superpose_with_options

__all__ = ["SuperposeError", "SuperposeOptions", "Superposition", "rmsd", "rmsd_flat", "superpose", "superpose_with_options"]

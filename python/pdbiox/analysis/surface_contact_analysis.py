"""Surface-contact analysis.

Rust names this module ``surface_contacts``. Python retains
``analysis.surface_contacts`` as the canonical callable, so the module surface
lives at this unambiguous path.
"""

from . import SurfaceContactOptions, surface_contacts

__all__ = ["SurfaceContactOptions", "surface_contacts"]

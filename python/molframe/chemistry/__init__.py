"""Chemical element data and per-atom radii."""
from .._native import chemistry as _native

ElementProperties = _native.ElementProperties
element = _native.element
vdw_radii = _native.vdw_radii
vdw_radius = _native.vdw_radius

__all__ = ["ElementProperties", "element", "vdw_radii", "vdw_radius"]

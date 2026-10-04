"""Numeric geometry kernels."""

from .._native import geometry as _native

Superposition = _native.Superposition
centroid = _native.centroid
distance_matrix = _native.distance_matrix
rmsd = _native.rmsd
distances = _native.distances
angles = _native.angles
dihedrals = _native.dihedrals
centre_of_mass = _native.centre_of_mass
radius_of_gyration = _native.radius_of_gyration
inertia_tensor = _native.inertia_tensor
principal_axes = _native.principal_axes
gyration_axes = _native.gyration_axes
asphericity = _native.asphericity
shape_parameter = _native.shape_parameter
best_fit_plane = _native.best_fit_plane
plane_deviation = _native.plane_deviation
rmsd_after_fit = _native.rmsd_after_fit
superpose = _native.superpose
rmsf = _native.rmsf

__all__ = [
    "Superposition",
    "angles",
    "asphericity",
    "best_fit_plane",
    "centre_of_mass",
    "centroid",
    "dihedrals",
    "distance_matrix",
    "distances",
    "gyration_axes",
    "inertia_tensor",
    "plane_deviation",
    "principal_axes",
    "radius_of_gyration",
    "rmsd",
    "rmsd_after_fit",
    "rmsf",
    "shape_parameter",
    "superpose",
]

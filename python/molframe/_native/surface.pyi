"""The compiled ``molframe._native.surface`` submodule."""

from .._io_types import (
    AtomDepthError, BuriedSurfaceError, SasaError, SurfaceGeometryError, SurfaceWorkflowError,
)
from .._trajectory_format_models import (
    SurfaceMesh,
)
from ..surface import (
    AtomContactArea, AtomDepthOptions, BuriedSurface, BuriedSurfaceOp, Cavity, CurvatureQuality,
    ExcludedSurfacePoint, IndexedSurfaceMesh, MeshReport, MoleculeRole, Sasa, SolventExcludedSurface,
    SurfaceComponent, SurfaceComponentError, SurfaceComponentFilter, SurfaceCurvature, SurfaceDistances, SurfaceFace,
    SurfaceGridOptions, SurfacePoint, SurfaceTriangle, atom_contact_areas, atom_depths, buried_solvent_excluded_surface,
    buried_solvent_excluded_surface_with_options, buried_surface, cavities, cavities_with_options, collect_shrake_rupley, edge_geodesic_distances,
    fibonacci_sphere, filter_surface_components, governed_surface_geometry, lee_richards, shrake_rupley, solvent_accessible_surface,
    solvent_excluded_surface, solvent_excluded_surface_with_options, surface_components, surface_curvatures, surface_patch, surface_points,
    surface_points_at_density, surface_points_excluding_pairs, visit_shrake_rupley, write_obj,
)
from ..surface._workflow import (
    SurfaceWorkflowOptions, SurfaceWorkflowResult, analyse_surface_geometry,
)

__all__: list[str]

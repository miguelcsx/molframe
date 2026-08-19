from .._native.traj import *
from .._native.traj import __all__
from .._native.traj import TrajectoryInterpolation, TrajectoryInterpolationError
from ._dimensionality import *
from .._trajectory import MeanSquaredDisplacementOp, PipelineReader, RmsdToReference
from .dielectric import DielectricEstimate, DielectricOptions, dielectric_from_dipoles
from . import format_xyz
from . import format_aims, format_charmm, format_gamess, format_gro, format_txyz
from . import format_amber
from . import format_dlpoly
from . import format_dcd, format_trr, format_xtc
from . import format_amber_netcdf
from . import format_gromacs_itp
from . import format_amber_topology, format_psf
from . import format_namd
from . import format_gromos, format_gsd, format_h5md, format_hoomd_xml, format_lammps, format_lammps_data, format_tng, format_tpr, format_trz
from . import (
    aims, amber, amber_netcdf, amber_topology, analysis,
    centroid_clustering, charmm, cluster_algorithms, clustering,
    contract_workflows, dcd, dielectric, diffusion, dimensionality,
    dihedral_pca, dispatch, displacement_statistics, dlpoly, dms,
    dynamic_query, electrostatics, ensemble, ensemble_similarity,
    ensemble_statistics, execution, format_dms, format_trc, frame_ops,
    frame_store, gamess, geometry_ensemble, governed, gro, gromacs_itp,
    gsd, h5md, hoomd_xml, imd, io_dispatch, lammps, lammps_data,
    manifold, minimal, msd, namd, neighbors, pca,
    periodic_ops, periodic_transform, protocol_imd, psf, reader, selection,
    similarity_ensemble, similarity_path, solvent_dynamics, spatial_frames,
    statistics_ensemble, tng, topology_amber, topology_minimal, tpr,
    trajectory, trajectory_stream, transform, trc, trr, trz, txyz,
    xtc, xyz,
)

def interpolate_trajectory_frames(
    frames: list[list[tuple[float, float, float]]],
    fraction: float,
    interpolation: TrajectoryInterpolation,
) -> list[tuple[float, float, float]]: ...

# Module/function collisions use centroid_clustering, similarity_path, and
# solvent_dynamics so the package-level native functions remain callable.

"""The compiled ``molframe._native.traj`` submodule."""

from numpy import float64
from numpy.typing import NDArray
from .._io_types import (
    DielectricError, DmsError, EnsembleGeometryError, EnsembleSimilarityError, EnsembleStatisticsError, GovernedEnsembleError,
    ImdError, KMeansError, MsdError, PathSimilarityError, TrajectoryError, TrajectoryIoError,
    UpdatingSelectionError, WaterDynamicsError,
)
from .._provider import (
    Frame,
)
from .._trajectory import (
    DEFAULT_PAIRWISE_MEMORY_LIMIT, DEFAULT_PROCRUSTES_MAXIMUM_ITERATIONS, DEFAULT_PROCRUSTES_TOLERANCE, CartesianFit,
)
from .._trajectory_dispatch import (
    AmberAsciiReadOptions, StreamFrame, TrajectoryReadOptions, TrajectoryReaderOptions, TrajectoryStreamReader, contact_counts_stream,
    read_trajectory, read_trajectory_materialized, rmsd_stream, rmsf_stream, run_analysis_stream, write_trajectory,
)
from .._trajectory_dms import (
    DmsBond, DmsCell, DmsFrame, DmsParticle, DmsSystem, DmsTopology,
    DmsVersion, read_dms, write_dms,
)
from .._trajectory_ensemble import (
    Clustering, ConvergenceBlock, EnsembleDistanceMatrix, FrameAlignment, GroupVariance, HarmonicSimilarity,
    HarmonicSimilarityOptions, KMeans, KMeansOptions, Linkage, RemainderPolicy, analyse_agglomerative_clustering,
    analyse_block_convergence, analyse_cluster_population_similarity, analyse_dbscan_clustering, analyse_generalized_procrustes_mean,
    analyse_group_coordinate_variance, analyse_harmonic_ensemble_similarity, analyse_kmeans, analyse_medoid,
    analyse_pairwise_fitted_rmsd, analyse_pairwise_torus_distance, analyse_rmsd_to_reference,
)
from .._trajectory_format_models import (
    AimsAtom, AimsError, AimsGeometry, AmberError, AmberRestart, AmberRestartLayout,
    CharmmAtom, CharmmCard, CharmmCardFormat, CharmmError, DiffusionMap, DlPolyAtom,
    DlPolyConfig, DlPolyError, DlPolyFrame, DlPolyHistory, FrameValue, GamessAtom,
    GamessError, GamessFrame, GamessRunType, GamessTrajectory, GroAtom, GroError,
    GroFrame, ImdClient, ImdConnectionOptions, ImdEnergies, ImdForce, ImdLimits,
    ImdMessage, ImdPeerEndian, PcaResult, PeriodicAngle, Rotation3, SurfaceMesh,
    Timestep, TxyzAtom, TxyzError,
    TxyzFrame, XyzAtom, XyzFrame, parse_aims_geometry, parse_amber_ascii_trajectory, parse_amber_restart,
    parse_amber_restart_record, parse_charmm_record, parse_dlpoly_config, parse_dlpoly_history, parse_gamess_output, parse_gro_records,
    parse_txyz_records, parse_xyz, surface_mesh, write_aims_geometry, write_amber_restart, write_charmm_card, write_dlpoly_config,
    write_dlpoly_history, write_gro, write_txyz, write_xyz,
)
from .._trajectory_kernels import (
    DielectricEstimate, MeanSquaredDisplacement, TrajectoryDielectricOptions, generalized_procrustes_mean, kmeans, mean_squared_displacement,
    pairwise_fitted_rmsd, pairwise_torus_distance, rmsd_to_reference,
)
from .._trajectory_operations import (
    GeneralizedProcrustesMean, MeanSquaredDisplacementOp, MeanSquaredDisplacementSeries, PairwiseFittedRmsd, RmsdToReference,
)
from .._trajectory_runtime import (
    Center, ChainedReader, Fit, FormatMetadata, FrameAnalysis, FrameNeighborList,
    MemoryReader, MinimalTopology, NeighborStatistics, PathFrameMetric, PathSimilarity, PipelineReader,
    RandomAccess, RigidTransform, StreamingReader, SurvivalMode, Trajectory, TrajectoryData,
    TrajectoryFormat, TrajectoryMetadata, TrajectoryUnits, TrajectoryWriteOptions, TrzWriteOptions, Units,
    Unwrap, UpdatingSelection, WaterDynamics, WaterSurvival, Wrap, agglomerative_clustering,
    analyse_cartesian_pca, analyse_diffusion_map, analyse_dihedral_pca, block_convergence, cartesian_pca, cluster_population_similarity,
    dbscan_clustering, diffusion_map, group_coordinate_variance, harmonic_ensemble_similarity, medoid, path_similarity,
    run_analysis, torsion_pca, trajectory_water_dynamics,
)
from ..analysis._analysis_operations import (
    water_dynamics,
)
from ..surface import (
    SurfaceGridOptions,
)
from ..traj import (
    interpolate_trajectory_frames,
)
from ..traj._dimensionality import (
    analyse_diffusion, analyse_pca, analyse_torsion_pca,
)
from ..traj.format_amber_netcdf import (
    AmberNetcdfError, AmberNetcdfMetadata, AmberNetcdfPrecision, AmberNetcdfTrajectory, AmberNetcdfWriteOptions, parse_amber_netcdf,
    parse_amber_netcdf_record, write_amber_netcdf,
)
from ..traj.format_amber_topology import (
    AmberSection, AmberTopology, AmberTopologyAngle, AmberTopologyAtom, AmberTopologyBond, AmberTopologyDihedral,
    AmberTopologyError, AmberTopologyResidue, parse_amber_topology,
)
from ..traj.format_dcd import (
    DcdEndian, DcdError, DcdHeader, DcdTrajectory, DcdWriteOptions, parse_dcd,
    write_dcd,
)
from ..traj.format_gromacs_itp import (
    GromacsInteraction, GromacsItp, GromacsItpAtom, GromacsItpError, GromacsMoleculeType, parse_gromacs_itp,
)
from ..traj.format_gromos import (
    GromosBoundary, GromosError, GromosTrajectory, parse_gromos11_trc,
)
from ..traj.format_gsd import (
    GsdError, GsdOptions, GsdTrajectory, parse_gsd, write_gsd,
)
from ..traj.format_h5md import (
    H5mdError, H5mdMetadata, H5mdOptions, H5mdTrajectory, H5mdUnitSystem, parse_h5md,
    parse_h5md_record_with_options, parse_h5md_with_options, write_h5md, write_h5md_with_metadata,
)
from ..traj.format_hoomd_xml import (
    HoomdBox, HoomdConfiguration, HoomdInteraction, HoomdXmlError, parse_hoomd_xml,
)
from ..traj.format_lammps import (
    LammpsError, parse_lammps_dump,
)
from ..traj.format_lammps_data import (
    LammpsAtomStyle, LammpsData, LammpsDataAtom, LammpsDataCell, LammpsDataError, LammpsInteraction,
    parse_lammps_data,
)
from ..traj.format_namd import (
    NamdBinary, NamdEndian, NamdError, parse_namd_binary, write_namd_binary,
)
from ..traj.format_psf import (
    PsfAtom, PsfError, PsfTopology, parse_psf, write_psf,
)
from ..traj.format_tng import (
    TngCompression, TngError, TngTrajectory, TngWriteOptions, parse_tng, write_tng,
)
from ..traj.format_tpr import (
    TprAtom, TprBond, TprError, TprHeader, TprResidue, TprTopology,
    parse_tpr,
)
from ..traj.format_trr import (
    TrrError, TrrPrecision, TrrTrajectory, TrrWriteOptions, parse_trr, write_trr,
    write_trr_with_precisions,
)
from ..traj.format_trz import (
    TrzError, TrzTrajectory, parse_trz, write_trz,
)
from ..traj.format_xtc import (
    XtcError, XtcTrajectory, XtcWriteOptions, parse_xtc, write_xtc, write_xtc_with_precisions,
)
from ..traj.xvg import (
    XvgSeries, read_xvg,
)

class TrajectoryInterpolation:
    Linear: TrajectoryInterpolation
    CentripetalCatmullRom: TrajectoryInterpolation

class TrajectoryInterpolationError:
    TopologyMismatch: TrajectoryInterpolationError
    NonFiniteInput: TrajectoryInterpolationError
    FractionOutOfRange: TrajectoryInterpolationError


def trajectory_dielectric_from_dipoles(
    dipoles: NDArray[float64], options: TrajectoryDielectricOptions
) -> DielectricEstimate: ...

__all__: list[str]

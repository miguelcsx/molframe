//! Trajectory frames and the formats that carry them.
//!
//! A trajectory is a sequence of coordinate frames over a fixed set of atoms.
//! This crate holds the in-memory frame model — random access to any frame and
//! slicing to a subset — and the readers that fill it. The coordinate store is
//! deliberately separate from any topology: frames are positions, and what those
//! positions belong to is decided a layer up.

#![forbid(unsafe_code)]

#[path = "kmeans.rs"]
pub mod centroid_clustering;
#[path = "clustering.rs"]
pub mod cluster_algorithms;
#[path = "governed/mod.rs"]
pub mod contract_workflows;
#[path = "pca.rs"]
pub mod dimensionality;
#[path = "msd.rs"]
pub mod displacement_statistics;
#[path = "selection.rs"]
pub mod dynamic_query;
#[path = "dielectric.rs"]
pub mod electrostatics;
#[path = "dcd_write.rs"]
mod encoder_dcd;
#[path = "dlpoly_write.rs"]
mod encoder_dlpoly;
#[path = "trr_write.rs"]
mod encoder_trr;
#[path = "trz_write.rs"]
mod encoder_trz;
#[path = "analysis.rs"]
pub mod execution;
#[path = "aims.rs"]
pub mod format_aims;
#[path = "amber.rs"]
pub mod format_amber;
#[path = "amber_netcdf/mod.rs"]
pub mod format_amber_netcdf;
#[path = "amber_restart/mod.rs"]
mod format_amber_restart;
#[path = "charmm.rs"]
pub mod format_charmm;
#[path = "dcd.rs"]
pub mod format_dcd;
#[path = "dlpoly.rs"]
pub mod format_dlpoly;
#[path = "dms/mod.rs"]
pub mod format_dms;
#[path = "gamess.rs"]
pub mod format_gamess;
#[path = "gro.rs"]
pub mod format_gro;
#[path = "gromacs_itp.rs"]
pub mod format_gromacs_itp;
#[path = "gsd.rs"]
pub mod format_gsd;
#[path = "h5md/mod.rs"]
pub mod format_h5md;
#[path = "hoomd_xml.rs"]
pub mod format_hoomd_xml;
#[path = "lammps.rs"]
pub mod format_lammps;
#[path = "lammps_data.rs"]
pub mod format_lammps_data;
#[path = "namd.rs"]
pub mod format_namd;
#[path = "psf.rs"]
pub mod format_psf;
#[path = "tng/mod.rs"]
pub mod format_tng;
#[path = "tpr/mod.rs"]
pub mod format_tpr;
#[path = "trc.rs"]
pub mod format_trc;
#[path = "trr.rs"]
pub mod format_trr;
#[path = "trz.rs"]
pub mod format_trz;
#[path = "txyz.rs"]
pub mod format_txyz;
#[path = "xtc.rs"]
pub mod format_xtc;
#[path = "xyz.rs"]
pub mod format_xyz;
#[path = "transform.rs"]
pub mod frame_ops;
#[path = "trajectory.rs"]
pub mod frame_store;
#[path = "frame_view.rs"]
mod frame_view;
#[path = "ensemble.rs"]
pub mod geometry_ensemble;
mod interpolation;
#[path = "dispatch/mod.rs"]
pub mod io_dispatch;
#[path = "diffusion.rs"]
pub mod manifold;
mod numeric;
#[path = "cell.rs"]
mod periodic_box;
#[path = "periodic_transform.rs"]
pub mod periodic_ops;
#[path = "imd/mod.rs"]
pub mod protocol_imd;
#[path = "ensemble_similarity.rs"]
pub mod similarity_ensemble;
#[path = "path_similarity.rs"]
pub mod similarity_path;
#[path = "water_dynamics.rs"]
pub mod solvent_dynamics;
#[path = "neighbors.rs"]
pub mod spatial_frames;
#[path = "ensemble_statistics.rs"]
pub mod statistics_ensemble;
#[path = "amber_topology.rs"]
pub mod topology_amber;
#[path = "minimal.rs"]
pub mod topology_minimal;
#[path = "reader.rs"]
pub mod trajectory_stream;

pub use centroid_clustering as kmeans;
pub use cluster_algorithms as clustering;
pub use contract_workflows as governed;
pub use dimensionality as pca;
pub use displacement_statistics as msd;
pub use dynamic_query as selection;
pub use electrostatics as dielectric;
use encoder_dcd as dcd_write;
use encoder_dlpoly as dlpoly_write;
use encoder_trr as trr_write;
use encoder_trz as trz_write;
pub use execution as analysis;
pub use format_aims as aims;
pub use format_amber as amber;
pub use format_amber_netcdf as amber_netcdf;
use format_amber_restart as amber_restart;
pub use format_charmm as charmm;
pub use format_dcd as dcd;
pub use format_dlpoly as dlpoly;
pub use format_dms as dms;
pub use format_gamess as gamess;
pub use format_gro as gro;
pub use format_gromacs_itp as gromacs_itp;
pub use format_gsd as gsd;
pub use format_h5md as h5md;
pub use format_hoomd_xml as hoomd_xml;
pub use format_lammps as lammps;
pub use format_lammps_data as lammps_data;
pub use format_namd as namd;
pub use format_psf as psf;
pub use format_tng as tng;
pub use format_tpr as tpr;
pub use format_trc as trc;
pub use format_trr as trr;
pub use format_trz as trz;
pub use format_txyz as txyz;
pub use format_xtc as xtc;
pub use format_xyz as xyz;
pub use frame_ops as transform;
pub use frame_store as trajectory;
pub use geometry_ensemble as ensemble;
pub use io_dispatch as dispatch;
pub use manifold as diffusion;
use periodic_box as cell;
pub use periodic_ops as periodic_transform;
pub use protocol_imd as imd;
pub use similarity_ensemble as ensemble_similarity;
pub use similarity_path as path_similarity;
pub use solvent_dynamics as water_dynamics;
pub use spatial_frames as neighbors;
pub use statistics_ensemble as ensemble_statistics;
pub use topology_amber as amber_topology;
pub use topology_minimal as minimal;
pub use trajectory_stream as reader;

pub use aims::{AimsAtom, AimsError, AimsGeometry, parse_aims_geometry, write_aims_geometry};
pub use amber::{
    AmberError, AmberRestartLayout, parse_amber_ascii_trajectory, parse_amber_restart,
    parse_amber_restart_record,
};
pub use amber_netcdf::{
    AmberNetcdfError, AmberNetcdfMetadata, AmberNetcdfPrecision, AmberNetcdfTrajectory,
    AmberNetcdfWriteOptions, parse_amber_netcdf, parse_amber_netcdf_record, write_amber_netcdf,
};
pub use amber_restart::{AmberRestart, write_amber_restart};
pub use amber_topology::{
    AmberSection, AmberTopology, AmberTopologyAngle, AmberTopologyAtom, AmberTopologyBond,
    AmberTopologyDihedral, AmberTopologyError, AmberTopologyResidue, parse_amber_topology,
};
pub use analysis::{FrameAnalysis, run_analysis};
pub use charmm::{
    CharmmAtom, CharmmCard, CharmmCardFormat, CharmmError, parse_charmm_record, write_charmm_card,
};
pub use clustering::{Clustering, Linkage, agglomerative_clustering, dbscan_clustering, medoid};
pub use dcd::{DcdEndian, DcdError, DcdHeader, DcdTrajectory, parse_dcd};
pub use dcd_write::{DcdWriteOptions, write_dcd};
pub use dielectric::{
    DielectricError, DielectricEstimate, DielectricOptions, dielectric_from_dipoles,
};
pub use diffusion::{DiffusionMap, diffusion_map};
pub use dispatch::{
    AmberAsciiReadOptions, FormatMetadata, TrajectoryData, TrajectoryFormat, TrajectoryIoError,
    TrajectoryMetadata, TrajectoryReadOptions, TrajectoryWriteOptions, TrzWriteOptions,
    read_trajectory, write_trajectory,
};
pub use dlpoly::{
    DlPolyAtom, DlPolyConfig, DlPolyError, DlPolyFrame, DlPolyHistory, parse_dlpoly_config,
    parse_dlpoly_history,
};
pub use dlpoly_write::{write_dlpoly_config, write_dlpoly_history};
pub use dms::{
    DmsBond, DmsCell, DmsError, DmsFrame, DmsParticle, DmsSystem, DmsTopology, DmsVersion,
    read_dms, write_dms,
};
pub use ensemble::{
    DEFAULT_PAIRWISE_MEMORY_LIMIT, DEFAULT_PROCRUSTES_MAXIMUM_ITERATIONS,
    DEFAULT_PROCRUSTES_TOLERANCE, EnsembleDistanceMatrix, EnsembleGeometryError, FrameAlignment,
    generalized_procrustes_mean, generalized_procrustes_mean_view, pairwise_fitted_rmsd,
    pairwise_fitted_rmsd_view, pairwise_torus_distance, rmsd_to_reference, rmsd_to_reference_view,
};
pub use ensemble_similarity::{
    EnsembleSimilarityError, HarmonicSimilarity, HarmonicSimilarityOptions,
    cluster_population_similarity, harmonic_ensemble_similarity,
};
pub use ensemble_statistics::{
    ConvergenceBlock, EnsembleStatisticsError, GroupVariance, RemainderPolicy, block_convergence,
    group_coordinate_variance, group_coordinate_variance_view,
};
pub use frame_view::{FrameView, FrameViewError};
pub use gamess::{
    GamessAtom, GamessError, GamessFrame, GamessRunType, GamessTrajectory, parse_gamess_output,
};
pub use governed::{
    GovernedEnsembleError, analyse_agglomerative_clustering, analyse_block_convergence,
    analyse_cartesian_pca, analyse_cartesian_pca_view, analyse_cluster_population_similarity,
    analyse_dbscan_clustering, analyse_diffusion_map, analyse_dihedral_pca,
    analyse_generalized_procrustes_mean, analyse_generalized_procrustes_mean_view,
    analyse_group_coordinate_variance, analyse_group_coordinate_variance_view,
    analyse_harmonic_ensemble_similarity, analyse_kmeans, analyse_kmeans_view,
    analyse_mean_squared_displacement_view, analyse_medoid, analyse_pairwise_fitted_rmsd,
    analyse_pairwise_fitted_rmsd_view, analyse_pairwise_torus_distance, analyse_rmsd_to_reference,
    analyse_rmsd_to_reference_view,
};
pub use gro::{GroAtom, GroError, GroFrame, parse_gro_records, write_gro};
pub use gromacs_itp::{
    GromacsInteraction, GromacsItp, GromacsItpAtom, GromacsItpError, GromacsMoleculeType,
    parse_gromacs_itp,
};
pub use gsd::{GsdError, GsdOptions, GsdTrajectory, parse_gsd, write_gsd};
pub use h5md::{
    H5mdError, H5mdMetadata, H5mdOptions, H5mdTrajectory, H5mdUnitSystem, parse_h5md,
    parse_h5md_record_with_options, parse_h5md_with_options, write_h5md, write_h5md_with_metadata,
};
pub use hoomd_xml::{
    HoomdBox, HoomdConfiguration, HoomdInteraction, HoomdXmlError, parse_hoomd_xml,
};
pub use imd::{
    ImdClient, ImdConnectionOptions, ImdEnergies, ImdError, ImdForce, ImdLimits, ImdMessage,
    ImdPeerEndian,
};
pub use interpolation::{
    TrajectoryInterpolation, TrajectoryInterpolationError, interpolate_trajectory_frames,
};
pub use kmeans::{KMeans, KMeansError, KMeansOptions, kmeans, kmeans_view};
pub use lammps::{LammpsError, parse_lammps_dump};
pub use lammps_data::{
    LammpsAtomStyle, LammpsData, LammpsDataAtom, LammpsDataCell, LammpsDataError,
    LammpsInteraction, parse_lammps_data,
};
pub use minimal::MinimalTopology;
pub use msd::{
    MeanSquaredDisplacement, MsdError, mean_squared_displacement, mean_squared_displacement_view,
};
pub use namd::{NamdBinary, NamdEndian, NamdError, parse_namd_binary, write_namd_binary};
pub use neighbors::{FrameNeighborList, NeighborStatistics};
pub use path_similarity::{
    PathFrameMetric, PathSimilarity, PathSimilarityError, path_similarity, path_similarity_view,
};
pub use pca::{CartesianFit, PcaResult, cartesian_pca, cartesian_pca_view, dihedral_pca};
pub use periodic_transform::{Unwrap, Wrap};
pub use psf::{PsfAtom, PsfError, PsfTopology, parse_psf, write_psf};
pub use reader::{
    ChainedReader, FrameValue, MemoryReader, RandomAccess, StreamingReader, Timestep,
    TrajectoryError, TrajectoryReader, Units,
};
pub use selection::{UpdatingSelection, UpdatingSelectionError};
pub use tng::{TngCompression, TngError, TngTrajectory, TngWriteOptions, parse_tng, write_tng};
pub use tpr::{TprAtom, TprBond, TprError, TprHeader, TprResidue, TprTopology, parse_tpr};
pub use trajectory::{Frame, Trajectory};
pub use transform::{Center, Fit, FrameTransform, PipelineReader, RigidTransform};
pub use trc::{GromosBoundary, GromosError, GromosTrajectory, parse_gromos11_trc};
pub use trr::{TrrError, TrrPrecision, TrrTrajectory, parse_trr};
pub use trr_write::{TrrWriteOptions, write_trr, write_trr_with_precisions};
pub use trz::{TrzError, TrzTrajectory, parse_trz};
pub use trz_write::write_trz;
pub use txyz::{TxyzAtom, TxyzError, TxyzFrame, parse_txyz_records, write_txyz};
pub use water_dynamics::{
    SurvivalMode, WaterDynamics, WaterDynamicsError, WaterSurvival, water_dynamics,
};
pub use xtc::{
    XtcError, XtcTrajectory, XtcWriteOptions, parse_xtc, write_xtc, write_xtc_with_precisions,
};
pub use xyz::{XyzAtom, XyzFrame, parse_xyz, write_xyz};

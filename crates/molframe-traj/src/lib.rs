//! Trajectory frames and the formats that carry them.
//!
//! A trajectory is a sequence of coordinate frames over a fixed set of atoms.
//! This crate holds the in-memory frame model — random access to any frame and
//! slicing to a subset — and the readers that fill it. The coordinate store is
//! deliberately separate from any topology: frames are positions, and what those
//! positions belong to is decided a layer up.

#![forbid(unsafe_code)]

pub mod analysis;
mod cell;
mod diagnostics;
pub mod dispatch;
mod ensembles;
mod formats;
mod frame_index;
mod frame_view;
pub mod imd;
mod interpolation;
mod numeric;
pub mod periodic_transform;
pub mod selection;
mod streaming;
pub mod trajectory;
pub mod transform;
pub mod trr;
mod words;

#[cfg(test)]
mod test_support;

pub use ensembles::{
    clustering, dielectric, diffusion, ensemble, ensemble_similarity, ensemble_statistics,
    governed, kmeans, msd, neighbors, path_similarity, pca, water_dynamics,
};
use ensembles::{pca_fit, pca_linalg, pca_memory};
pub use formats::{
    aims, amber, amber_netcdf, amber_topology, charmm, dcd, dlpoly, dms, gamess, gro, gromacs_itp,
    gsd, h5md, hoomd_xml, lammps, lammps_data, minimal, namd, psf, tng, tpr, trc, trz, txyz, xtc,
    xvg, xyz,
};
use formats::{amber_restart, dcd_write, dlpoly_write, trr_write, trz_write};
pub use streaming::reader;
use streaming::{
    batch, ordered_blocks, selected, stream_chunks, stream_consume, stream_contacts,
    stream_metrics, stream_pipeline, stream_ranges,
};

// Former module names, kept so every existing path still resolves.
pub use analysis as execution;
pub use dispatch as io_dispatch;
pub use ensembles::clustering as cluster_algorithms;
pub use ensembles::dielectric as electrostatics;
pub use ensembles::diffusion as manifold;
pub use ensembles::ensemble as geometry_ensemble;
pub use ensembles::ensemble_similarity as similarity_ensemble;
pub use ensembles::ensemble_statistics as statistics_ensemble;
pub use ensembles::governed as contract_workflows;
pub use ensembles::kmeans as centroid_clustering;
pub use ensembles::msd as displacement_statistics;
pub use ensembles::neighbors as spatial_frames;
pub use ensembles::path_similarity as similarity_path;
pub use ensembles::pca as dimensionality;
pub use ensembles::water_dynamics as solvent_dynamics;
pub use formats::aims as format_aims;
pub use formats::amber as format_amber;
pub use formats::amber_netcdf as format_amber_netcdf;
pub use formats::amber_topology as topology_amber;
pub use formats::charmm as format_charmm;
pub use formats::dcd as format_dcd;
pub use formats::dlpoly as format_dlpoly;
pub use formats::dms as format_dms;
pub use formats::gamess as format_gamess;
pub use formats::gro as format_gro;
pub use formats::gromacs_itp as format_gromacs_itp;
pub use formats::gsd as format_gsd;
pub use formats::h5md as format_h5md;
pub use formats::hoomd_xml as format_hoomd_xml;
pub use formats::lammps as format_lammps;
pub use formats::lammps_data as format_lammps_data;
pub use formats::minimal as topology_minimal;
pub use formats::namd as format_namd;
pub use formats::psf as format_psf;
pub use formats::tng as format_tng;
pub use formats::tpr as format_tpr;
pub use formats::trc as format_trc;
pub use formats::trz as format_trz;
pub use formats::txyz as format_txyz;
pub use formats::xtc as format_xtc;
pub use formats::xvg as format_xvg;
pub use formats::xyz as format_xyz;
pub use imd as protocol_imd;
pub use periodic_transform as periodic_ops;
pub use selection as dynamic_query;
pub use streaming::reader as trajectory_stream;
pub use trajectory as frame_store;
pub use transform as frame_ops;
pub use trr as format_trr;

pub use format_trr::TrrReader;
pub use selected::{DenseChunk, FrameSlice, SelectedReader, count_frames, read_dense_chunk};
pub use stream_chunks::ParallelChunks;
pub use stream_consume::run_analysis_stream;
pub use stream_contacts::contact_counts_stream;
pub use stream_metrics::{rmsd_stream, rmsf_stream};
pub use stream_pipeline::run_frame_pipeline;
pub use stream_ranges::run_frame_ranges;

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
pub use batch::{TrajectoryBatch, TrajectoryBatchSource};
pub use charmm::{
    CharmmAtom, CharmmCard, CharmmCardFormat, CharmmError, parse_charmm_record, write_charmm_card,
};
pub use clustering::{Clustering, Linkage, agglomerative_clustering, dbscan_clustering, medoid};
pub use dcd::{DcdEndian, DcdError, DcdHeader, DcdReader, DcdTrajectory, parse_dcd};
pub use dcd_write::{DcdWriteOptions, write_dcd};
pub use dielectric::{
    DielectricError, DielectricEstimate, DielectricOptions, dielectric_from_dipoles,
};
pub use diffusion::{DiffusionMap, diffusion_map};
pub use dispatch::{
    AmberAsciiReadOptions, DensePositions, FormatMetadata, TrajectoryData, TrajectoryFormat,
    TrajectoryIoError, TrajectoryMetadata, TrajectoryReadOptions, TrajectoryReaderOptions,
    TrajectoryWriteOptions, TrzWriteOptions, read_trajectory, read_trajectory_in,
    read_trajectory_materialized, write_trajectory,
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
pub use frame_index::{FileStamp, FrameIndex, FrameIndexError};
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
pub use periodic_transform::{TemporalUnwrap, Unwrap, Wrap};
pub use psf::{PsfAtom, PsfError, PsfTopology, parse_psf, write_psf};
pub use reader::{
    ChainedReader, FrameValue, MemoryReader, RandomAccess, StreamingReader, Timestep,
    TrajectoryError, TrajectoryReader, Units,
};
pub use selection::{UpdatingSelection, UpdatingSelectionError};
pub use tng::{TngCompression, TngError, TngTrajectory, TngWriteOptions, parse_tng, write_tng};
pub use tpr::{TprAtom, TprBond, TprError, TprHeader, TprResidue, TprTopology, parse_tpr};
pub use trajectory::{Frame, FrameRef, Trajectory, TrajectoryBuildError};
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
pub use words::UnknownWord;
pub use xtc::{
    XtcError, XtcReader, XtcTrajectory, XtcWriteOptions, parse_xtc, write_xtc,
    write_xtc_with_precisions,
};
pub use xyz::{XyzAtom, XyzFrame, parse_xyz, write_xyz};

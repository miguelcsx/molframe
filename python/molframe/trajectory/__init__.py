"""Trajectory operations."""

from .._native import trajectory as _native

Clustering = _native.Clustering
Convergence = _native.Convergence
DiffusionMap = _native.DiffusionMap
GroupVariance = _native.GroupVariance
HarmonicSimilarity = _native.HarmonicSimilarity
KMeans = _native.KMeans
MeanSquaredDisplacement = _native.MeanSquaredDisplacement
PathSimilarity = _native.PathSimilarity
Pca = _native.Pca
Trajectory = _native.Trajectory
agglomerative_clustering = _native.agglomerative_clustering
block_convergence = _native.block_convergence
dbscan = _native.dbscan
diffusion_map = _native.diffusion_map
group_variance = _native.group_variance
harmonic_similarity = _native.harmonic_similarity
kmeans = _native.kmeans
mean_structure = _native.mean_structure
medoid = _native.medoid
msd = _native.msd
pairwise_rmsd = _native.pairwise_rmsd
path_similarity = _native.path_similarity
pca = _native.pca
population_similarity = _native.population_similarity
read = _native.read
rmsd = _native.rmsd
FrameSeries = _native.FrameSeries
TrajectoryChunks = _native.TrajectoryChunks
contact_counts_stream = _native.contact_counts_stream
sasa_stream = _native.sasa_stream
read_chunks = _native.read_chunks
rmsd_stream = _native.rmsd_stream
rmsf_stream = _native.rmsf_stream
run_analysis_stream = _native.run_analysis_stream

__all__ = [
    "Clustering",
    "Convergence",
    "DiffusionMap",
    "FrameSeries",
    "GroupVariance",
    "HarmonicSimilarity",
    "KMeans",
    "MeanSquaredDisplacement",
    "PathSimilarity",
    "Pca",
    "Trajectory",
    "TrajectoryChunks",
    "agglomerative_clustering",
    "block_convergence",
    "contact_counts_stream",
    "dbscan",
    "diffusion_map",
    "group_variance",
    "harmonic_similarity",
    "kmeans",
    "mean_structure",
    "medoid",
    "msd",
    "pairwise_rmsd",
    "path_similarity",
    "pca",
    "population_similarity",
    "read",
    "read_chunks",
    "rmsd",
    "rmsd_stream",
    "rmsf_stream",
    "run_analysis_stream",
    "sasa_stream",
]

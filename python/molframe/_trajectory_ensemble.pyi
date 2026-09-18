from typing import Any, final
from numpy import float32, float64
from numpy.typing import NDArray
from .core.contract import Analysis
from .query import AnalysisPolicy

@final
class FrameAlignment:
    Unaligned: FrameAlignment
    Rigid: FrameAlignment

@final
class Linkage:
    Single: Linkage
    Complete: Linkage
    Average: Linkage

@final
class RemainderPolicy:
    Include: RemainderPolicy
    Reject: RemainderPolicy

@final
class EnsembleDistanceMatrix:
    def __init__(self, values: NDArray[float64]) -> None: ...
    size: int
    values: list[float]

@final
class Clustering:
    labels: list[int | None]
    members: list[list[int]]
    medoids: list[int]

@final
class KMeansOptions:
    def __init__(
        self,
        initial_centres: list[int],
        maximum_iterations: int,
        convergence_tolerance_squared: float,
    ) -> None: ...

@final
class KMeans:
    labels: list[int]
    centres: list[list[float]]
    inertia: float
    iterations: int

@final
class HarmonicSimilarityOptions:
    def __init__(self, covariance_regularization: float, memory_limit_bytes: int) -> None: ...

@final
class HarmonicSimilarity:
    mean_term: float
    covariance_term: float
    distance: float
    similarity: float

@final
class GroupVariance:
    atoms: list[int]
    variance_by_axis: tuple[float, float, float]
    rms_fluctuation: float

@final
class ConvergenceBlock:
    start: int
    end: int
    block_mean: float
    cumulative_mean: float

def analyse_rmsd_to_reference(frames: NDArray[float32], reference: int, alignment: FrameAlignment, policy: AnalysisPolicy) -> Analysis[Any]: ...
def analyse_pairwise_fitted_rmsd(frames: NDArray[float32], memory_limit: int, policy: AnalysisPolicy) -> Analysis[Any]: ...
def analyse_pairwise_torus_distance(points: NDArray[float64], weights: NDArray[float64], metric_name: str, memory_limit: int, policy: AnalysisPolicy) -> Analysis[Any]: ...
def analyse_generalized_procrustes_mean(frames: NDArray[float32], tolerance: float, maximum_iterations: int, policy: AnalysisPolicy) -> Analysis[Any]: ...
def analyse_agglomerative_clustering(distances: EnsembleDistanceMatrix, cluster_count: int, linkage: Linkage, policy: AnalysisPolicy) -> Analysis[Any]: ...
def analyse_dbscan_clustering(distances: EnsembleDistanceMatrix, epsilon: float, minimum_points: int, policy: AnalysisPolicy) -> Analysis[Any]: ...
def analyse_medoid(distances: EnsembleDistanceMatrix, members: list[int], policy: AnalysisPolicy) -> Analysis[Any]: ...
def analyse_kmeans(observations: NDArray[float64], options: KMeansOptions, policy: AnalysisPolicy) -> Analysis[Any]: ...
def analyse_harmonic_ensemble_similarity(first: list[list[float]], second: list[list[float]], options: HarmonicSimilarityOptions, policy: AnalysisPolicy) -> Analysis[Any]: ...
def analyse_cluster_population_similarity(first: list[int], second: list[int], cluster_count: int, policy: AnalysisPolicy) -> Analysis[Any]: ...
def analyse_group_coordinate_variance(frames: NDArray[float32], groups: list[list[int]], policy: AnalysisPolicy) -> Analysis[Any]: ...
def analyse_block_convergence(values: list[float], block_size: int, remainder: RemainderPolicy, policy: AnalysisPolicy) -> Analysis[Any]: ...

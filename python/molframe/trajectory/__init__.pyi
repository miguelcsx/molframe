from collections.abc import Sequence
from os import PathLike
from typing import Literal, Protocol, final

from .. import Analysis, AnalysisPolicy, Float32Array, Float64Array

class Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

@final
class Trajectory:
    @property
    def format(self) -> str: ...
    @property
    def positions(self) -> Array: ...
    @property
    def times(self) -> Array: ...
    @property
    def n_frames(self) -> int: ...
    @property
    def n_atoms(self) -> int: ...
    def __len__(self) -> int: ...

def read(
    path: str | PathLike[str],
    *,
    format: Literal["xtc", "trr", "dcd", "tng", "gro", "xyz", "lammps_dump", "netcdf"]
    | None = None,
) -> Trajectory: ...
def rmsd(
    positions: Array,
    *,
    reference: int = 0,
    align: bool = True,
    policy: AnalysisPolicy | None = None,
) -> Analysis[Array]: ...

class UInt64Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class Int64Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

@final
class Pca:
    """Principal components of a set of observations."""

    @property
    def mean(self) -> Float64Array:
        """The mean of every feature."""
    @property
    def eigenvalues(self) -> Float64Array:
        """Eigenvalues, largest first."""
    @property
    def components(self) -> Float64Array:
        """One axis per component, ``(components, features)``, each with a canonical sign."""
    @property
    def projections(self) -> Float64Array:
        """Each observation's coordinates on the components, ``(observations, components)``."""

@final
class DiffusionMap:
    """A diffusion-map embedding of an affinity geometry."""

    @property
    def eigenvalues(self) -> Float64Array: ...
    @property
    def coordinates(self) -> Float64Array:
        """Observation coordinates after the diffusion time, ``(observations, dimensions)``."""
    @property
    def graph_components(self) -> UInt64Array:
        """The connected component of the affinity graph each observation lies in."""

@final
class KMeans:
    """A k-means partition."""

    @property
    def labels(self) -> UInt64Array: ...
    @property
    def centres(self) -> Float64Array: ...
    @property
    def inertia(self) -> float: ...
    @property
    def iterations(self) -> int: ...

@final
class Clustering:
    """A partition into clusters, each with a medoid representative."""

    def __len__(self) -> int: ...
    @property
    def labels(self) -> Int64Array:
        """The cluster of each observation; ``-1`` marks DBSCAN noise."""
    @property
    def members(self) -> list[list[int]]: ...
    @property
    def medoids(self) -> list[int]: ...

@final
class MeanSquaredDisplacement:
    """Mean squared displacement by frame lag."""

    def __len__(self) -> int: ...
    @property
    def lag(self) -> UInt64Array: ...
    @property
    def observations(self) -> UInt64Array: ...
    @property
    def value(self) -> Float64Array: ...

@final
class Convergence:
    """Mean per-block and cumulative values of a series."""

    def __len__(self) -> int: ...
    @property
    def start(self) -> UInt64Array: ...
    @property
    def end(self) -> UInt64Array: ...
    @property
    def block_mean(self) -> Float64Array: ...
    @property
    def cumulative_mean(self) -> Float64Array: ...

@final
class GroupVariance:
    """The coordinate variance of one atom group."""

    @property
    def atoms(self) -> list[int]: ...
    @property
    def variance_by_axis(self) -> tuple[float, float, float]: ...
    @property
    def rms_fluctuation(self) -> float: ...

@final
class HarmonicSimilarity:
    """A distance between two Gaussian models of an ensemble."""

    @property
    def mean_term(self) -> float: ...
    @property
    def covariance_term(self) -> float: ...
    @property
    def distance(self) -> float: ...
    @property
    def similarity(self) -> float:
        """``exp(-distance)``: one for identical models."""

@final
class PathSimilarity:
    @property
    def hausdorff_distance(self) -> float: ...
    @property
    def discrete_frechet_distance(self) -> float: ...

def pairwise_rmsd(
    positions: Array,
    *,
    memory_limit: int = 1073741824,
    policy: AnalysisPolicy | None = None,
) -> Analysis[Float64Array]:
    """Return the RMSD between every two frames after the best rigid fit of each pair.

    ``positions`` has shape ``(frames, atoms, 3)``; the matrix is ``(frames, frames)``.
    """

def pca(
    positions: Array,
    *,
    components: int,
    fit: Literal["none", "mean"] = "none",
    reference: Array | None = None,
    tolerance: float = 1e-6,
    max_iterations: int = 100,
    memory_limit: int = 1073741824,
    policy: AnalysisPolicy | None = None,
) -> Analysis[Pca]:
    """Return the principal components of the Cartesian coordinates of every frame.

    ``fit`` aligns frames first: ``"none"`` keeps them, ``"mean"`` fits to an iterative
    Procrustes mean, and a ``reference`` of shape ``(atoms, 3)`` fits to that structure.
    The covariance is normalised by the number of frames less one.
    """

def mean_structure(
    positions: Array,
    *,
    tolerance: float = 1e-6,
    max_iterations: int = 100,
    policy: AnalysisPolicy | None = None,
) -> Analysis[Float32Array]:
    """Return the mean structure of the frames from generalized Procrustes alignment."""

def diffusion_map(
    distances: Float64Array,
    *,
    epsilon: float,
    dimensions: int,
    time: int = 1,
    policy: AnalysisPolicy | None = None,
) -> Analysis[DiffusionMap]:
    """Return a diffusion map of a square distance matrix, affinities ``exp(-d^2 / epsilon)``."""

def kmeans(
    observations: Float64Array,
    *,
    initial_centres: Sequence[int],
    maximum_iterations: int = 300,
    tolerance_squared: float = 1e-12,
    policy: AnalysisPolicy | None = None,
) -> Analysis[KMeans]:
    """Return k-means on ``(observations, features)``, started from the listed observations."""

def agglomerative_clustering(
    distances: Float64Array,
    *,
    clusters: int,
    linkage: Literal["single", "complete", "average"] = "average",
    policy: AnalysisPolicy | None = None,
) -> Analysis[Clustering]:
    """Return an agglomerative clustering of a square distance matrix."""

def dbscan(
    distances: Float64Array,
    *,
    epsilon: float,
    minimum_points: int,
    policy: AnalysisPolicy | None = None,
) -> Analysis[Clustering]:
    """Return DBSCAN on a square distance matrix; the rest are noise (label ``-1``)."""

def medoid(
    distances: Float64Array,
    members: Sequence[int],
    *,
    policy: AnalysisPolicy | None = None,
) -> Analysis[int]:
    """Return the member with the least total distance to the other members."""

def msd(
    positions: Array,
    *,
    maximum_lag: int,
    atoms: Sequence[int] | None = None,
    policy: AnalysisPolicy | None = None,
) -> Analysis[MeanSquaredDisplacement]:
    """Return window-averaged mean squared displacement for lags ``0..=maximum_lag``.

    Coordinates are used as given: unwrap a periodic trajectory first.
    """

def group_variance(
    positions: Array,
    groups: Sequence[Sequence[int]],
    *,
    policy: AnalysisPolicy | None = None,
) -> Analysis[list[GroupVariance]]:
    """Return the coordinate variance of each atom group over frames already aligned.

    The variance of a group is the mean over its atoms of each atom's population variance
    about its own mean position.
    """

def block_convergence(
    values: Sequence[float],
    *,
    block_size: int,
    remainder: Literal["include", "reject"] = "include",
    policy: AnalysisPolicy | None = None,
) -> Analysis[Convergence]:
    """Return block means and cumulative means of a series."""

def harmonic_similarity(
    first: Float64Array,
    second: Float64Array,
    *,
    regularization: float,
    memory_limit: int = 1073741824,
    policy: AnalysisPolicy | None = None,
) -> Analysis[HarmonicSimilarity]:
    """Return the distance between Gaussian models fitted to two sets of observations.

    Each model's covariance is normalised by its observation count, and ``regularization``
    is added to every diagonal element.
    """

def population_similarity(
    first: Sequence[int],
    second: Sequence[int],
    *,
    clusters: int,
    policy: AnalysisPolicy | None = None,
) -> Analysis[float]:
    """Return the Jensen-Shannon similarity of two ensembles' cluster populations."""

def path_similarity(
    first: Array,
    second: Array,
    *,
    metric: Literal["fitted_rmsd", "cartesian_rmsd"] = "fitted_rmsd",
    memory_limit: int = 1073741824,
) -> PathSimilarity:
    """Return the Hausdorff and discrete Frechet distances between two paths of frames."""

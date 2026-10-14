from collections.abc import Callable, Sequence
from os import PathLike
from typing import Literal, Protocol, final

from .. import (
    Analysis,
    AnalysisPolicy,
    ExecutionContext,
    Float32Array,
    Float64Array,
    Structure,
)
from ..crystal import UnitCell

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
    def cells(self) -> Array | None:
        """`(frames, 6)` rows of a, b, c, alpha, beta, gamma; NaN where a frame has no cell."""
    @property
    def frame_indices(self) -> UInt64Array:
        """Source frame of each row; sliced and strided reads keep the original numbering."""
    def cell(self, frame: int) -> UnitCell | None: ...
    @property
    def n_frames(self) -> int: ...
    @property
    def n_atoms(self) -> int: ...
    def __len__(self) -> int: ...

def read(
    path: str | PathLike[str],
    *,
    format: Literal[
        "xtc",
        "trr",
        "dcd",
        "tng",
        "gro",
        "xyz",
        "lammps_dump",
        "netcdf",
        "pdb",
        "mmcif",
        "bcif",
    ]
    | None = None,
    frames: slice | None = None,
    atoms: Sequence[int] | None = None,
    memory_limit: int | None = None,
) -> Trajectory:
    """Read a trajectory, or a slice and atom subset of it, into one exactly sized buffer.

    A ``.pdb``, ``.cif``, ``.mmcif`` or ``.bcif`` path reads a multi-model
    structure and returns its models as frames. A ragged ensemble has no shared
    atom axis and raises an error naming ``ragged_models``.

    `frames` is a positive-stride `slice` of source frames and `atoms` the atom
    indices to keep. XTC seeks; DCD and TRR decode and discard in between.
    """

@final
class TrajectoryChunks:
    """Blocks of frames from one trajectory; memory stays at one block."""

    def __iter__(self) -> TrajectoryChunks: ...
    def __next__(self) -> Trajectory: ...
    def seek(self, frame: int) -> None:
        """Move to the `frame`-th selected frame, forward or back."""
    @property
    def n_frames(self) -> int | None: ...
    @property
    def n_atoms(self) -> int: ...
    @property
    def chunk_frames(self) -> int: ...

def read_chunks(
    path: str | PathLike[str],
    *,
    chunk_frames: int = 256,
    format: Literal["xtc", "trr", "dcd", "tng", "gro", "xyz", "lammps_dump", "netcdf"]
    | None = None,
    frames: slice | None = None,
    atoms: Sequence[int] | None = None,
    memory_limit: int | None = None,
    context: ExecutionContext | None = None,
) -> TrajectoryChunks:
    """Iterate blocks of frames.

    `context=` with more than one worker decodes ahead in parallel, in file
    order, within the context's memory budget.
    """

@final
class FrameSeries:
    """One value per streamed frame, with its source frame and time."""

    @property
    def frames(self) -> UInt64Array: ...
    @property
    def times(self) -> Float64Array: ...
    @property
    def values(self) -> Float64Array | UInt64Array: ...
    def __len__(self) -> int: ...

def rmsd_stream(
    path: str | PathLike[str],
    reference: int | Array,
    *,
    align: bool = True,
    format: Literal["xtc", "trr", "dcd", "tng", "gro", "xyz", "lammps_dump", "netcdf"]
    | None = None,
    frames: slice | None = None,
    atoms: Sequence[int] | None = None,
    memory_limit: int | None = None,
    context: ExecutionContext | None = None,
) -> FrameSeries:
    """RMSD of every frame against a reference array or source frame, in constant memory."""

def rmsf_stream(
    path: str | PathLike[str],
    *,
    format: Literal["xtc", "trr", "dcd", "tng", "gro", "xyz", "lammps_dump", "netcdf"]
    | None = None,
    frames: slice | None = None,
    atoms: Sequence[int] | None = None,
    memory_limit: int | None = None,
    context: ExecutionContext | None = None,
) -> Float64Array:
    """Per-atom RMSF in constant memory; no alignment is applied."""

def contact_counts_stream(
    path: str | PathLike[str],
    cutoff: float,
    *,
    backend: str = "auto",
    format: Literal["xtc", "trr", "dcd", "tng", "gro", "xyz", "lammps_dump", "netcdf"]
    | None = None,
    frames: slice | None = None,
    atoms: Sequence[int] | None = None,
    memory_limit: int | None = None,
    context: ExecutionContext | None = None,
) -> FrameSeries:
    """Count unique atom pairs within `cutoff` in every frame, in constant memory."""

def run_analysis_stream(
    path: str | PathLike[str],
    callback: Callable[[int, float | None, Float32Array], object],
    *,
    format: Literal["xtc", "trr", "dcd", "tng", "gro", "xyz", "lammps_dump", "netcdf"]
    | None = None,
    frames: slice | None = None,
    atoms: Sequence[int] | None = None,
    memory_limit: int | None = None,
    context: ExecutionContext | None = None,
) -> int:
    """Call `callback(frame, time, positions)` per frame; the GIL is free between calls."""

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

def sasa_stream(
    path: str | PathLike[str],
    topology: Structure,
    *,
    probe: float | None = None,
    method: str | None = None,
    n_points: int | None = None,
    hetero: str | None = None,
    metals: str | None = None,
    radii: str | None = None,
    policy: AnalysisPolicy | None = None,
    format: Literal["xtc", "trr", "dcd", "tng", "gro", "xyz", "lammps_dump", "netcdf"]
    | None = None,
    frames: slice | None = None,
    memory_limit: int | None = None,
    context: ExecutionContext | None = None,
) -> Analysis[FrameSeries]:
    """Total solvent-accessible area of every frame, in square angstroms, in constant memory.

    ``topology`` supplies the atoms and their radii, so the trajectory must have its atoms in
    the same order (a different count raises ``E5102``). The options, the radius rule and the
    policy are those of :func:`molframe.analysis.sasa`, including ``E7103`` for a ``radii``
    that conflicts with the policy; only the Shrake-Rupley method streams. The values are the
    per-atom areas of the batch calculation summed over the same atoms; hydrogens, other
    conformers and atoms with no radius are neither measured nor read from the trajectory.
    """

"""Analysis result and option types shared across native analysis kernels."""

from typing import Any, final
from numpy import bool_, float32, float64
from numpy.typing import NDArray
from ..geom import EigenOptions
from .. import Structure
from ..chem import ComponentDictionary, StereoConfiguration
from ..core import MissingResidue
from ..core.contract import Analysis
from .._io_types import (BasePairError, CationPiError, DensityError, DsspBinaryError, DsspError,
    DynamicsError, EnsembleGeometryError, EnsembleSimilarityError,
    EnsembleStatisticsError, FragmentMappingError, GnmError, GovernedAnalysisError,
    GovernedEnsembleError, GovernedNativeError, HelicalError, HseError,
    HydrogenBondError, KMeansError, NativeError, NucleicTorsionError,
    PhysicalKernelError, PiStackingError, PolymerError, PoreError, RadialError,
    StandaloneAnalysisError)
from ..query import AnalysisPolicy, Namespace, Selection

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
    def __init__(self, initial_centres: list[int], maximum_iterations: int, convergence_tolerance_squared: float) -> None: ...

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

@final
class Contact:
    def __init__(self, first: int, second: int, distance: float) -> None: ...
    first: int
    second: int
    distance: float
    def __repr__(self) -> str: ...

@final
class BasePairOptions:
    def __init__(self, maximum_donor_acceptor_distance: float, minimum_angle_degrees: float, minimum_hydrogen_bonds: int, backend: object, periodic: bool) -> None: ...

@final
class BasePair:
    first: int
    second: int
    hydrogen_bond_count: int
    closest_distance: float

@final
class ResidueContact:
    first: int
    second: int
    min_distance: float

@final
class ContactMap:
    residue_count: int
    contacts: list[ResidueContact]

@final
class FragmentReference:
    def __init__(self, id: str, coordinates_array: NDArray[float32]) -> None: ...

@final
class FragmentMatch:
    start: int
    fragment_id: str
    rmsd: float

@final
class GnmOptions:
    def __init__(self, contact_distance: float, mode_count: int, zero_mode_tolerance: float, memory_limit_bytes: int, backend: object, reduction: object = ...) -> None: ...
    @property
    def reduction(self) -> object: ...

@final
class GaussianNetworkModel:
    sites: NDArray
    eigenvalues: NDArray[float64]
    modes: NDArray[float64]
    zero_modes: int

@final
class AnmOptions:
    def __init__(self, contact_distance: float, mode_count: int, zero_mode_tolerance: float, memory_limit_bytes: int, backend: object) -> None: ...

@final
class AnisotropicNetworkModel:
    sites: NDArray
    eigenvalues: NDArray[float64]
    modes: NDArray[float64]
    zero_modes: int
    def fluctuations(self) -> NDArray[float64]: ...
    def project(self, values: NDArray[float64]) -> NDArray[float64]: ...
    def displace(self, positions: NDArray[float32], amplitudes: NDArray[float64]) -> NDArray[float32]: ...

@final
class NormalMode:
    def __init__(self, index: int, scale: float, displacements: NDArray[float64]) -> None: ...
    index: int
    scale: float
    displacements: NDArray[float64]

@final
class NormalModeSet:
    name: str | None
    atom_names: list[str]
    residue_names: list[str]
    residue_ids: list[int]
    chain_ids: list[str]
    b_factors: list[float]
    coordinates: NDArray[float64]
    modes: list[NormalMode]
    def __len__(self) -> int: ...
    def displace(self, amplitudes: NDArray[float64]) -> NDArray[float64]: ...

@final
class HalfSphereExposure:
    residue: int
    upper: int
    lower: int

@final
class NativeContacts:
    native: int
    kept: int
    fraction: float

@final
class NucleicTorsions:
    residue: int
    alpha: float | None
    beta: float | None
    gamma: float | None
    delta: float | None
    epsilon: float | None
    zeta: float | None
    chi: float | None

@final
class Pucker:
    phase_degrees: float
    amplitude: float

@final
class SurfaceContactOptions:
    def __init__(self, tolerance: float, probe: float, surface_density: float, minimum_area: float, backend: object) -> None: ...

@final
class SseKind:
    AlphaHelix: SseKind
    Strand: SseKind
    Turn: SseKind
    Coil: SseKind

@final
class SecondaryStructure:
    residue: int
    kind: SseKind

@final
class DsspOptions:
    def __init__(self, electrostatic_prefactor: float, hydrogen_bond_energy: float, amide_hydrogen_distance: float, minimum_sequence_separation: int, helix_offset: int, turn_offsets: tuple[int, int]) -> None: ...

@final
class DsspSegment:
    kind: str
    begin_chain: str
    begin_sequence: int
    end_chain: str
    end_sequence: int

@final
class LinearDensityOptions:
    def __init__(self, axis: CartesianAxis, minimum: float, maximum: float, bins: int) -> None: ...
    axis: CartesianAxis
    minimum: float
    maximum: float
    bins: int

@final
class SseRecord:
    residue: int
    kind: SseKind

@final
class CentreGroup:
    def __init__(self, atoms: Selection) -> None: ...
    atoms: Selection

@final
class LeafletOptions:
    def __init__(self, connection_distance: float, *, backend: object = ...) -> None: ...

@final
class Leaflet:
    sites: list[int]

@final
class QualityIssue:
    ZeroOccupancy: QualityIssue
    OccupancyOutOfRange: QualityIssue
    NegativeBFactor: QualityIssue

@final
class QualityFlag:
    atom: int
    issue: QualityIssue

@final
class BondDeviation:
    atom_a: int
    atom_b: int
    observed: float
    expected: float
    deviation: float

@final
class Clash:
    first: int
    second: int
    overlap: float

@final
class ChainCompleteness:
    chain: str
    observed: int
    canonical: int
    missing: list[MissingResidue]

@final
class CisPeptide:
    residue: int
    omega: float

@final
class PlanarityFlag:
    residue: int
    deviation: float

@final
class PlanarityOptions:
    def __init__(self, maximum_deviation: float, plane_fit: object) -> None: ...

@final
class ValenceError:
    atom: int
    bonds: int
    maximum: int

@final
class ReferenceAssessment:
    set: str
    version: str
    distribution: str
    probability: float
    percentile: float | None

@final
class RamachandranRegion:
    AlphaHelixRight: RamachandranRegion
    BetaSheet: RamachandranRegion
    AlphaHelixLeft: RamachandranRegion
    Outlier: RamachandranRegion

@final
class RamachandranRecord:
    residue: int
    phi: float
    psi: float
    region: RamachandranRegion
    reference_id: str
    reference_version: str
    distribution: str
    probability: float
    percentile: float | None

@final
class ReferenceDistribution:
    @staticmethod
    def histogram(name: str, edges: list[float], weights: list[float]) -> ReferenceDistribution: ...
    @staticmethod
    def grid(name: str, x_edges: list[float], y_edges: list[float], weights: list[float]) -> ReferenceDistribution: ...

@final
class ReferenceLibrary:
    def __init__(self, id: str, version: str, distributions: list[ReferenceDistribution]) -> None: ...
    id: str
    version: str

@final
class RamachandranBasin:
    @staticmethod
    def alpha_right(distribution: str) -> RamachandranBasin: ...
    @staticmethod
    def beta_sheet(distribution: str) -> RamachandranBasin: ...
    @staticmethod
    def alpha_left(distribution: str) -> RamachandranBasin: ...

@final
class RamachandranOptions:
    def __init__(self, references: ReferenceLibrary, basins: list[RamachandranBasin], minimum_probability: float) -> None: ...

@final
class ChiralityIssue:
    Inverted: ChiralityIssue
    Degenerate: ChiralityIssue

@final
class ChiralityOptions:
    def __init__(self, minimum_abs_volume: float) -> None: ...

@final
class ChiralityFlag:
    residue: int
    centre: int
    expected: StereoConfiguration
    issue: ChiralityIssue
    observed_volume: float
    reference_volume: float

@final
class ChiralityReport:
    flags: list[ChiralityFlag]
    findings: list[str]
    dictionary_version: str

@final
class RotamerDefinition:
    def __init__(self, component: str, chi_index: int, atoms: list[str], distribution: str) -> None: ...

@final
class RotamerProfile:
    def __init__(self, id: str, version: str, definitions: list[RotamerDefinition]) -> None: ...
    id: str
    version: str

@final
class RotamerOptions:
    def __init__(self, minimum_probability: float) -> None: ...

@final
class RotamerFlag:
    residue: int
    chi_index: int
    chi_degrees: float
    probability: float
    distribution: str

@final
class RotamerReport:
    flags: list[RotamerFlag]
    findings: list[str]
    dictionary_version: str
    profile_id: str
    profile_version: str
    reference_id: str
    reference_version: str

@final
class PolymerStatistics:
    contour_length: float
    end_to_end_distance: float
    mean_segment_length: float
    adjacent_tangent_correlation: float | None
    persistence_length: float | None

@final
class PoreOptions:
    def __init__(self, axis: tuple[list[float], list[float]], start: float, end: float, samples: int, search_radius: float, grid_spacing: float, probe_radius: float, *, memory_limit_bytes: int = ...) -> None: ...

PoreProfileOptions = PoreOptions

@final
class PoreSample:
    axial_coordinate: float
    centre: list[float]
    radius: float

@final
class RadialOptions:
    def __init__(self, minimum: float, maximum: float, bins: int, volume: float, *, backend: object = ...) -> None: ...

RadialDistributionOptions = RadialOptions

@final
class RadialBin:
    lower: float
    upper: float
    count: int
    distribution: float

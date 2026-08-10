from typing import final
from numpy import bool_, float32, float64
from numpy.typing import NDArray
from .geometry import EigenOptions
from . import Structure
from .chemistry import ComponentDictionary
from .contract import Analysis
from .query import AnalysisPolicy, Namespace, Selection
from .validation import (BFactorDistribution, BFactorOutlier, TlsBFactorFlag,
    TlsBFactorReport, TlsGroup, TlsModel, analyse_b_factor_distribution,
    analyse_tls_b_factor_consistency, b_factor_distribution, tls_b_factor_consistency)

@final
class Contact:
    first: int
    second: int
    distance: float

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
    def __init__(self, contact_distance: float, mode_count: int, zero_mode_tolerance: float, memory_limit_bytes: int, backend: object) -> None: ...

@final
class GaussianNetworkModel:
    sites: NDArray
    eigenvalues: NDArray[float64]
    modes: NDArray[float64]
    zero_modes: int

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
class MissingResidue:
    canonical_position: int
    component: str

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
class StereoConfiguration:
    R: StereoConfiguration
    S: StereoConfiguration
    Mixed: StereoConfiguration

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
    def __init__(self, axis: tuple[list[float], list[float]], start: float, end: float, samples: int, search_radius: float, grid_spacing: float, probe_radius: float) -> None: ...

@final
class PoreSample:
    axial_coordinate: float
    centre: list[float]
    radius: float

@final
class RadialOptions:
    def __init__(self, minimum: float, maximum: float, bins: int, volume: float, *, backend: object = ...) -> None: ...

@final
class RadialBin:
    lower: float
    upper: float
    count: int
    distribution: float

@final
class HydrogenBondOptions:
    def __init__(self, maximum_distance: float, minimum_angle: float, *, backend: object = ..., periodic: bool = False) -> None: ...

@final
class HydrogenBond:
    donor: int
    hydrogen: int
    acceptor: int
    donor_acceptor_distance: float
    hydrogen_acceptor_distance: float
    angle: float

@final
class SaltBridge:
    anion: int
    cation: int
    distance: float

@final
class StackingKind:
    Parallel: StackingKind
    TShaped: StackingKind

@final
class PiStacking:
    first: int
    second: int
    centre_distance: float
    angle: float
    kind: StackingKind

@final
class PiStackingOptions:
    def __init__(self, maximum_centre_distance: float, maximum_parallel_angle: float, minimum_perpendicular_angle: float, plane_fit: EigenOptions) -> None: ...

@final
class CationPi:
    cation_residue: int
    ring_residue: int
    distance: float

@final
class CationPiOptions:
    def __init__(self, maximum_distance: float, maximum_face_angle: float, plane_fit: EigenOptions) -> None: ...

@final
class WaterBridge:
    water: int
    first: int
    second: int

@final
class CartesianAxis:
    X: CartesianAxis
    Y: CartesianAxis
    Z: CartesianAxis

@final
class LinearDensityBin:
    lower: float
    upper: float
    weight: float
    density: float

@final
class DensityGridSpec:
    def __init__(self, origin: list[float], spacing: list[float], shape: list[int]) -> None: ...
    origin: list[float]
    spacing: list[float]
    shape: list[int]

@final
class DensityGrid:
    spec: DensityGridSpec
    density: NDArray[float64]
    excluded_weight: float

@final
class SurfaceAreas:
    first_alone: float
    second_alone: float
    together: float
    buried: float

@final
class AtomDepthOptions:
    def __init__(self, cell_size: float) -> None: ...

@final
class SurfaceGridOptions:
    STANDARD_MAX_CELLS: int
    def __init__(self, resolution: float, max_cells: int) -> None: ...
    @staticmethod
    def standard(resolution: float) -> SurfaceGridOptions: ...

@final
class Cavity:
    volume: float
    representative: list[float]
    cells: int

def polymer_statistics(path: NDArray[float32]) -> PolymerStatistics: ...
def pore_profile(positions: NDArray[float32], radii: NDArray[float32], options: PoreOptions) -> list[PoreSample]: ...
def linear_density(positions: NDArray[float32], weights: NDArray[float64], axis: CartesianAxis, bounds: tuple[float, float], bins: int) -> list[LinearDensityBin]: ...
def density_map(positions: NDArray[float32], weights: NDArray[float64], spec: DensityGridSpec) -> DensityGrid: ...
def solvent_accessible_surface(positions: NDArray[float32], radii: NDArray[float32], probe_radius: float, sample_points: int) -> list[float]: ...
def buried_surface(positions: NDArray[float32], radii: NDArray[float32], first: NDArray[bool_], probe_radius: float, sample_points: int) -> SurfaceAreas: ...
def atom_depths(atoms: NDArray[float32], surface: NDArray[float32], options: AtomDepthOptions) -> NDArray[float32]: ...
def cavities(positions: NDArray[float32], radii: NDArray[float32], probe: float, options: SurfaceGridOptions) -> list[Cavity]: ...
def assess_bond_deviation(deviation: BondDeviation, references: ReferenceLibrary, distribution: str) -> ReferenceAssessment: ...
def assess_ramachandran(record: RamachandranRecord, references: ReferenceLibrary, distribution: str) -> ReferenceAssessment: ...
def classify_ramachandran(phi: float, psi: float, options: RamachandranOptions) -> tuple[RamachandranRegion, ReferenceAssessment]: ...
def map_fragments(trace: list[tuple[float, float, float] | None], library: list[FragmentReference], maximum_rmsd: float) -> list[FragmentMatch]: ...
def sugar_pucker(torsions: tuple[float, float, float, float, float]) -> Pucker: ...

@final
class WaterDynamicsOptions:
    def __init__(self, maximum_lag: int) -> None: ...
    maximum_lag: int
@final
class WaterLag:
    lag: int; observations: int; surviving: int; resident: int
    survival_probability: float | None; residence_probability: float | None
@final
class DielectricOptions:
    def __init__(self, volume: float, temperature: float, fluctuation_prefactor: float) -> None: ...
    volume: float; temperature: float; fluctuation_prefactor: float
@final
class DielectricResult:
    mean_dipole: list[float]; fluctuation: float; relative_permittivity: float
@final
class BaseFrame:
    def __init__(self, origin: list[float], x: list[float], y: list[float], z: list[float]) -> None: ...
    origin: list[float]; x: list[float]; y: list[float]; z: list[float]
@final
class HelicalOptions:
    def __init__(self, frame_tolerance: float) -> None: ...
    frame_tolerance: float
@final
class HelicalParameters:
    x_displacement: float; y_displacement: float; z_displacement: float
    x_rotation_degrees: float; y_rotation_degrees: float; z_rotation_degrees: float
@final
class AltlocOccupancyOptions:
    def __init__(self, expected_sum: float, tolerance: float) -> None: ...
    expected_sum: float; tolerance: float
@final
class AltlocOccupancyIssue:
    MissingOccupancy: AltlocOccupancyIssue; SumMismatch: AltlocOccupancyIssue
@final
class AltlocOccupancyRecord:
    residue: int; atom_name: str; alternatives: int; assessed: int
    occupancy_sum: float | None; issue: AltlocOccupancyIssue | None
@final
class AltlocOccupancyReport:
    intended: int; assessed: int; records: list[AltlocOccupancyRecord]
@final
class ResidueAtomCompleteness:
    residue: int; component: str; intended: int; assessed: int
    missing: list[str]; ambiguous: list[str]
@final
class CcdCompletenessReport:
    intended: int; assessed: int; ambiguous: int; residues: list[ResidueAtomCompleteness]
@final
class PlaneRestraint:
    def __init__(self, id: str, atoms: Selection) -> None: ...
    id: str; atoms: Selection
@final
class PlaneRestraintFlag:
    id: str; deviation: float
@final
class PlaneRestraintReport:
    intended: int; assessed: int; flags: list[PlaneRestraintFlag]

def water_dynamics(occupancy: list[Selection], options: WaterDynamicsOptions) -> list[WaterLag]: ...
def analyse_water_dynamics(occupancy: list[Selection], options: WaterDynamicsOptions, policy: AnalysisPolicy) -> Analysis[list[WaterLag]]: ...
def dielectric_from_dipoles(dipoles: list[list[float]], options: DielectricOptions) -> DielectricResult: ...
def analyse_dielectric_from_dipoles(dipoles: list[list[float]], options: DielectricOptions, policy: AnalysisPolicy) -> Analysis[DielectricResult]: ...
def helical_parameters(first: BaseFrame, second: BaseFrame, options: HelicalOptions) -> HelicalParameters: ...
def helical_steps(frames: list[BaseFrame], options: HelicalOptions) -> list[HelicalParameters]: ...
def analyse_helical_parameters(first: BaseFrame, second: BaseFrame, options: HelicalOptions, policy: AnalysisPolicy) -> Analysis[HelicalParameters]: ...
def analyse_helical_steps(frames: list[BaseFrame], options: HelicalOptions, policy: AnalysisPolicy) -> Analysis[list[HelicalParameters]]: ...
def analyse_centre_of_mass_radial_distribution(structure: Structure, masses: list[float], left: list[Selection], right: list[Selection], options: RadialOptions, policy: AnalysisPolicy) -> Analysis[list[RadialBin]]: ...
def validate_altloc_occupancy(structure: Structure, namespace: Namespace, options: AltlocOccupancyOptions) -> AltlocOccupancyReport: ...
def analyse_altloc_occupancy(structure: Structure, options: AltlocOccupancyOptions, policy: AnalysisPolicy) -> Analysis[AltlocOccupancyReport]: ...
def validate_ccd_completeness(structure: Structure, dictionary: ComponentDictionary, policy: AnalysisPolicy) -> CcdCompletenessReport: ...
def analyse_ccd_completeness(structure: Structure, dictionary: ComponentDictionary, policy: AnalysisPolicy) -> Analysis[CcdCompletenessReport]: ...
def validate_plane_restraints(structure: Structure, restraints: list[PlaneRestraint], options: PlanarityOptions) -> PlaneRestraintReport: ...
def analyse_plane_restraints(structure: Structure, restraints: list[PlaneRestraint], options: PlanarityOptions, policy: AnalysisPolicy) -> Analysis[PlaneRestraintReport]: ...

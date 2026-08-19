from typing import Any, final
from numpy import bool_, float64
from numpy.typing import NDArray
from .._native.validate import *
from .._native.validate import __all__
from ._governed import *
from ._thermal_motion import (
    BFactorDistribution, BFactorOutlier, TlsBFactorFlag, TlsBFactorReport,
    TlsGroup, TlsModel,
)
from ..chem import ComponentDictionary, PolymerRoleProfile, RadiusSet, StereoConfiguration
from ..core import MissingResidue, Structure
from ..spatial import SpatialBackend
from ..core.contract import Analysis
from ..geom import EigenOptions
from ..query import AnalysisPolicy, Namespace, Selection
from ..xtal.density import DensityMap, MapBoundary
from .._io_types import (
    AltlocOccupancyError, AltlocOccupancyKernelError, BFactorError,
    BFactorKernelError, CcdCompletenessKernelError, CompletenessError,
    GovernedMapError, LigandGeometryKernelError, NucleicGeometryError,
    PlanarityError, PlaneRestraintKernelError, RamachandranError,
    RealSpaceCorrelationError, ReferenceError, RotamerError,
)

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
class AltlocOccupancyOptions:
    def __init__(self, expected_sum: float, tolerance: float) -> None: ...
    expected_sum: float
    tolerance: float

@final
class AltlocOccupancyIssue:
    MissingOccupancy: AltlocOccupancyIssue
    SumMismatch: AltlocOccupancyIssue

@final
class AltlocOccupancyRecord:
    residue: int
    atom_name: str
    alternatives: int
    assessed: int
    occupancy_sum: float | None
    issue: AltlocOccupancyIssue | None

@final
class AltlocOccupancyReport:
    intended: int
    assessed: int
    records: list[AltlocOccupancyRecord]

@final
class ResidueAtomCompleteness:
    residue: int
    component: str
    intended: int
    assessed: int
    missing: list[str]
    ambiguous: list[str]

@final
class CcdCompletenessReport:
    intended: int
    assessed: int
    ambiguous: int
    residues: list[ResidueAtomCompleteness]

@final
class PlaneRestraint:
    def __init__(self, id: str, atoms: Selection) -> None: ...
    id: str
    atoms: Selection

@final
class PlaneRestraintFlag:
    id: str
    deviation: float

@final
class PlaneRestraintReport:
    intended: int
    assessed: int
    flags: list[PlaneRestraintFlag]

class _ValidationOperation:
    def execute(self, structure: Structure) -> Analysis[Any]: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...

@final
class StericClashes(_ValidationOperation):
    def __init__(self, *, tolerance: float = ..., radii: RadiusSet | None = None, backend: SpatialBackend | None = None, policy: AnalysisPolicy | None = ...) -> None: ...
    def execute(self, structure: Structure) -> Analysis[Any]: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> StericClashes: ...
    def __repr__(self) -> str: ...

@final
class BondLengthDeviations(_ValidationOperation):
    def __init__(self, *, tolerance: float, policy: AnalysisPolicy | None = ...) -> None: ...
    def execute(self, structure: Structure) -> Analysis[Any]: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> BondLengthDeviations: ...
    def __repr__(self) -> str: ...

@final
class CisPeptides(_ValidationOperation):
    def __init__(self, *, threshold_degrees: float, policy: AnalysisPolicy | None = ...) -> None: ...
    def execute(self, structure: Structure) -> Analysis[Any]: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> CisPeptides: ...
    def __repr__(self) -> str: ...

@final
class PlanarityCheck(_ValidationOperation):
    def __init__(self, *, options: PlanarityOptions, policy: AnalysisPolicy | None = ...) -> None: ...
    def execute(self, structure: Structure) -> Analysis[Any]: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> PlanarityCheck: ...
    def __repr__(self) -> str: ...

@final
class QualityFlags(_ValidationOperation):
    def __init__(self, *, policy: AnalysisPolicy | None = ...) -> None: ...
    def execute(self, structure: Structure) -> Analysis[Any]: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> QualityFlags: ...
    def __repr__(self) -> str: ...

@final
class Valence(_ValidationOperation):
    def __init__(self, *, policy: AnalysisPolicy | None = ...) -> None: ...
    def execute(self, structure: Structure) -> Analysis[Any]: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> Valence: ...
    def __repr__(self) -> str: ...

@final
class Completeness(_ValidationOperation):
    def __init__(self, *, policy: AnalysisPolicy | None = ...) -> None: ...
    def execute(self, structure: Structure) -> Analysis[Any]: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> Completeness: ...
    def __repr__(self) -> str: ...

@final
class LigandGeometryReport:
    outliers: list[BondDeviation]
    intended: int
    assessed: int

@final
class RealSpaceCorrelation:
    coefficient: float
    sample_count: int
    observed_mean: float
    calculated_mean: float

@final
class ReferenceGeometryOptions:
    def __init__(self, maximum_bond_deviation: float, maximum_angle_deviation_degrees: float) -> None: ...
    maximum_bond_deviation: float
    maximum_angle_deviation_degrees: float

@final
class ReferenceBondFlag:
    residue: int
    first: int
    second: int
    observed: float
    expected: float
    deviation: float

@final
class ReferenceAngleFlag:
    residue: int
    first: int
    centre: int
    third: int
    observed_degrees: float
    expected_degrees: float
    deviation_degrees: float

@final
class ReferenceGeometryReport:
    bonds: list[ReferenceBondFlag]
    angles: list[ReferenceAngleFlag]
    findings: list[str]
    intended: int
    assessed: int
    options: ReferenceGeometryOptions

@final
class NucleicGeometryPolicy:
    def __init__(self, maximum_base_plane_deviation: float, glycosidic_bond_range: tuple[float, float], phosphodiester_bond_range: tuple[float, float], plane_fit: EigenOptions) -> None: ...
    maximum_base_plane_deviation: float
    glycosidic_bond_range: tuple[float, float]
    phosphodiester_bond_range: tuple[float, float]
    plane_fit: EigenOptions

@final
class NucleicGeometryIssue:
    kind: str
    available: int | None
    expected: int | None
    deviation: float | None
    distance: float | None

@final
class NucleicGeometryRecord:
    residue: int
    torsions: NucleicTorsions
    pucker: Pucker | None
    base_plane_deviation: float | None
    glycosidic_bond_length: float | None
    incoming_phosphodiester_length: float | None
    issues: list[NucleicGeometryIssue]

def altloc_occupancy_sums(structure: Structure, namespace: Namespace, options: AltlocOccupancyOptions) -> AltlocOccupancyReport: ...
def bond_length_deviations(structure: Structure, tolerance: float) -> list[BondDeviation]: ...
def ccd_missing_atoms(structure: Structure, dictionary: ComponentDictionary, policy: AnalysisPolicy) -> CcdCompletenessReport: ...
def chirality_outliers(structure: Structure, dictionary: ComponentDictionary, options: ChiralityOptions, *, policy: AnalysisPolicy | None = ...) -> ChiralityReport: ...
def cis_peptides(structure: Structure, threshold_degrees: float) -> list[CisPeptide]: ...
def clashes(structure: Structure, tolerance: float, radius_set: RadiusSet, *, backend: object = ...) -> list[Clash]: ...
def classify(phi: float, psi: float, options: RamachandranOptions) -> tuple[RamachandranRegion, ReferenceAssessment]: ...
def completeness(structure: Structure, namespace: Namespace) -> list[ChainCompleteness]: ...
def governed_masked_real_space_correlation(observed: DensityMap, calculated: DensityMap, mask: NDArray[bool_], policy: AnalysisPolicy) -> Analysis[RealSpaceCorrelation]: ...
def governed_real_space_map_correlation(observed: DensityMap, calculated: DensityMap, policy: AnalysisPolicy) -> Analysis[RealSpaceCorrelation]: ...
def governed_sampled_real_space_correlation(observed: DensityMap, calculated: DensityMap, positions: NDArray[float64], policy: AnalysisPolicy) -> Analysis[RealSpaceCorrelation]: ...
def ligand_geometry(structure: Structure, tolerance: float) -> LigandGeometryReport: ...
def ligand_geometry_outliers(structure: Structure, tolerance: float) -> list[BondDeviation]: ...
def masked_real_space_correlation(observed: DensityMap, calculated: DensityMap, mask: NDArray[bool_]) -> RealSpaceCorrelation: ...
def nonplanar_aromatic_rings(structure: Structure, options: PlanarityOptions) -> list[PlanarityFlag]: ...
def nucleic_acid_geometry(structure: Structure, provider: ComponentDictionary, roles: PolymerRoleProfile, policy: NucleicGeometryPolicy) -> list[NucleicGeometryRecord]: ...
def overvalent_atoms(structure: Structure) -> list[ValenceError]: ...
def plane_restraint_outliers(structure: Structure, restraints: list[PlaneRestraint], options: PlanarityOptions) -> PlaneRestraintReport: ...
def quality_flags(structure: Structure) -> list[QualityFlag]: ...
def ramachandran(structure: Structure, options: RamachandranOptions, *, outliers_only: bool = ...) -> list[RamachandranRecord]: ...
def ramachandran_outliers(structure: Structure, options: RamachandranOptions) -> list[RamachandranRecord]: ...
def real_space_map_correlation(observed: DensityMap, calculated: DensityMap) -> RealSpaceCorrelation: ...
def reference_geometry(structure: Structure, provider: ComponentDictionary, namespace: Namespace, options: ReferenceGeometryOptions) -> ReferenceGeometryReport: ...
def rotamer_outliers(structure: Structure, dictionary: ComponentDictionary, references: ReferenceLibrary, profile: RotamerProfile, options: RotamerOptions, *, policy: AnalysisPolicy | None = ...) -> RotamerReport: ...
def sampled_real_space_correlation(observed: DensityMap, calculated: DensityMap, positions: NDArray[float64], boundary: MapBoundary) -> RealSpaceCorrelation: ...
def validate_altloc_occupancy(structure: Structure, namespace: Namespace, options: AltlocOccupancyOptions) -> AltlocOccupancyReport: ...
def analyse_altloc_occupancy(structure: Structure, options: AltlocOccupancyOptions, policy: AnalysisPolicy) -> Analysis[AltlocOccupancyReport]: ...
def validate_ccd_completeness(structure: Structure, dictionary: ComponentDictionary, policy: AnalysisPolicy) -> CcdCompletenessReport: ...
def analyse_ccd_completeness(structure: Structure, dictionary: ComponentDictionary, policy: AnalysisPolicy) -> Analysis[CcdCompletenessReport]: ...
def validate_plane_restraints(structure: Structure, restraints: list[PlaneRestraint], options: PlanarityOptions) -> PlaneRestraintReport: ...
def analyse_plane_restraints(structure: Structure, restraints: list[PlaneRestraint], options: PlanarityOptions, policy: AnalysisPolicy) -> Analysis[PlaneRestraintReport]: ...
def analyse_b_factor_distribution(structure: Structure, selection: Selection, outlier_standard_deviations: float, policy: AnalysisPolicy) -> Analysis[BFactorDistribution]: ...
def analyse_tls_b_factor_consistency(structure: Structure, groups: list[TlsGroup], maximum_absolute_deviation: float, symmetry_tolerance: float, policy: AnalysisPolicy) -> Analysis[TlsBFactorReport]: ...

from os import PathLike
from typing import Mapping, Sequence
from numpy.typing import NDArray
from numpy import float32
from ..chem import ComponentDictionary
from ..geom import EigenOptions, Rigid
from ..query import AnalysisPolicy
from .. import Structure

class MotifError(ValueError): ...
class MappingError(ValueError): ...
class MeasurementError(ValueError): ...
class EvaluationError(ValueError): ...
class SpecificationError(ValueError): ...
class AmeMeasurementError(ValueError): ...

class AtomSite:
    def __init__(self, component: str, atom: str) -> None: ...
    component: str
    atom: str
    def __str__(self) -> str: ...

class ComponentRole:
    Residue: ComponentRole
    Ligand: ComponentRole
    Cofactor: ComponentRole
    Metal: ComponentRole
    Unknown: ComponentRole

class ComponentSpec:
    def __init__(self, role: ComponentRole) -> None: ...
    @staticmethod
    def residue() -> ComponentSpec: ...
    @staticmethod
    def ligand() -> ComponentSpec: ...
    @staticmethod
    def cofactor() -> ComponentSpec: ...
    @staticmethod
    def metal() -> ComponentSpec: ...
    def component(self, id: str) -> ComponentSpec: ...
    def require(self, atom: str) -> ComponentSpec: ...
    def equivalent(self, atoms: Sequence[str]) -> ComponentSpec: ...
    role: ComponentRole
    component_ids: list[str]
    required_atoms: list[str]
    equivalent_atoms: list[list[str]]

class Constraint:
    @staticmethod
    def distance(first: AtomSite, second: AtomSite, target: float, tolerance: float) -> Constraint: ...
    @staticmethod
    def angle(atoms: Sequence[AtomSite], target: float, tolerance: float) -> Constraint: ...
    @staticmethod
    def dihedral(atoms: Sequence[AtomSite], target: float, tolerance: float) -> Constraint: ...
    @staticmethod
    def chirality(atoms: Sequence[AtomSite], positive: bool) -> Constraint: ...
    @staticmethod
    def planarity(atoms: Sequence[AtomSite], tolerance: float) -> Constraint: ...
    @staticmethod
    def coordination(centre: AtomSite, partners: Sequence[AtomSite], count: int, max_distance: float) -> Constraint: ...
    @staticmethod
    def steric_exclusion(first: AtomSite, second: AtomSite, min_distance: float) -> Constraint: ...
    kind: str
    atoms: list[AtomSite]
    target: float | None
    tolerance: float | None
    positive: bool | None
    count: int | None
    max_distance: float | None
    min_distance: float | None

class NamedConstraint:
    def __init__(self, name: str, constraint: Constraint) -> None: ...
    name: str
    constraint: Constraint

class Motif:
    def __init__(self, components: Mapping[str, ComponentSpec], constraints: Sequence[NamedConstraint]) -> None: ...
    components: dict[str, ComponentSpec]
    constraints: list[NamedConstraint]
    def __repr__(self) -> str: ...

class AlignmentKind:
    NotRequired: AlignmentKind
    CallerSupplied: AlignmentKind

class MappedMotif:
    def __init__(self, components: Mapping[str, int], atoms: Sequence[tuple[AtomSite, Sequence[int]]]) -> None: ...
    components: dict[str, int]
    atoms: list[tuple[AtomSite, list[int]]]
    def atom_indices(self, site: AtomSite) -> list[int] | None: ...

class MappingSet:
    def __init__(self, mappings: Sequence[MappedMotif], ambiguous: bool) -> None: ...
    mappings: list[MappedMotif]
    ambiguous: bool

class AlignedMotif:
    def __init__(self, mapping: MappedMotif, transform: Rigid, kind: AlignmentKind) -> None: ...
    mapping: MappedMotif
    transform: Rigid
    kind: AlignmentKind

class Comparison:
    @staticmethod
    def less_than(limit: float) -> Comparison: ...
    @staticmethod
    def at_most(limit: float) -> Comparison: ...
    @staticmethod
    def greater_than(limit: float) -> Comparison: ...
    @staticmethod
    def at_least(limit: float) -> Comparison: ...
    @staticmethod
    def between(low: float, high: float) -> Comparison: ...
    kind: str
    values: list[float]

class MissingVerdict:
    Indeterminate: MissingVerdict
    Fail: MissingVerdict

class VerdictRule:
    def __init__(self, metric: str, comparison: Comparison) -> None: ...
    metric: str
    comparison: Comparison

class VerdictProfile:
    def __init__(self, id: str, rules: Sequence[VerdictRule], missing: MissingVerdict) -> None: ...
    id: str
    rules: list[VerdictRule]
    def decide(self, metrics: Mapping[str, float]) -> Verdict: ...

class RuleOutcome:
    metric: str
    value: float | None
    passed: bool | None

class VerdictStatus:
    Pass: VerdictStatus
    Fail: VerdictStatus
    Indeterminate: VerdictStatus

class Verdict:
    profile: str
    status: VerdictStatus
    outcomes: list[RuleOutcome]

class MeasurementValue:
    @staticmethod
    def scalar(value: float) -> MeasurementValue: ...
    @staticmethod
    def chirality(sign: int) -> MeasurementValue: ...
    @staticmethod
    def count(value: int) -> MeasurementValue: ...
    kind: str
    numeric: float
    scalar_value: float | None
    chirality_sign: int | None
    count_value: int | None

class IndeterminateReason:
    MissingCoordinate: IndeterminateReason
    DegenerateGeometry: IndeterminateReason

class ConstraintMeasurement:
    name: str
    value: MeasurementValue | None
    deviation: float | None
    satisfied: bool | None
    atoms: list[int]
    indeterminate: IndeterminateReason | None

class MeasurementSet:
    constraints: list[ConstraintMeasurement]
    metrics: dict[str, float]

class MeasurementOptions:
    def __init__(self, maximum_alternatives: int, *, plane_fit: EigenOptions | None = ...) -> None: ...
    maximum_alternatives: int
    plane_fit: EigenOptions

class Evaluation:
    mapping_index: int
    alignment: AlignmentKind
    measurements: MeasurementSet
    verdict: Verdict

# The name the extension registers this class under; the ``fx`` namespace binds
# ``Evaluation`` to it, so the two are the same object under two names.
FxEvaluation = Evaluation

class EvaluationReport:
    mapping_ambiguous: bool
    evaluations: list[Evaluation]

class ProfileAtomSet:
    FullScaffoldCAlpha: ProfileAtomSet
    MotifBackboneWithOxygen: ProfileAtomSet
    CatalyticBackbone: ProfileAtomSet
    CatalyticHeavyAtoms: ProfileAtomSet
    LigandAndBackbone: ProfileAtomSet
    Unknown: ProfileAtomSet

class ProfileAlignment:
    MeasuredAtoms: ProfileAlignment
    CatalyticBackbone: ProfileAlignment
    PredictionFrame: ProfileAlignment
    Unknown: ProfileAlignment

class ProfileMetric:
    name: str
    atoms: ProfileAtomSet
    alignment: ProfileAlignment

class CandidateAggregation:
    Any: CandidateAggregation

class CompatibilityProfile:
    id: str
    metrics: list[ProfileMetric]
    aggregation: CandidateAggregation
    def decide_candidate(self, metrics: Mapping[str, float]) -> Verdict: ...
    def decide_candidates(self, candidates: Sequence[Mapping[str, float]]) -> CompatibilityVerdict: ...

class CompatibilityVerdict:
    profile: str
    status: VerdictStatus
    candidates: list[Verdict]

class MotifBenchMeasurements:
    rmsd: float
    motif_rmsd: float
    def metrics(self) -> dict[str, float]: ...

class AmeMeasurements:
    catalytic_heavy_atom_rmsd: float
    ligand_backbone_min_distance: float
    def metrics(self) -> dict[str, float]: ...

class EvaluationSpecification:
    motif: Motif
    profile: VerdictProfile

def align_intrinsic(mapping: MappedMotif) -> AlignedMotif: ...
def align_with_transform(mapping: MappedMotif, transform: Rigid) -> AlignedMotif: ...
def map_motif(structure: Structure, motif: Motif, limit: int, *, dictionary: ComponentDictionary | None = ..., policy: AnalysisPolicy | None = ...) -> MappingSet: ...
def measure_constraints(structure: Structure, motif: Motif, aligned: AlignedMotif, options: MeasurementOptions) -> MeasurementSet: ...
def evaluate_motif(structure: Structure, motif: Motif, profile: VerdictProfile, mapping_limit: int, measurement: MeasurementOptions, *, dictionary: ComponentDictionary | None = ..., policy: AnalysisPolicy | None = ...) -> EvaluationReport: ...
def measure_motifbench(generated_c_alpha: NDArray[float32], predicted_c_alpha: NDArray[float32], reference_motif_backbone: NDArray[float32], predicted_motif_backbone: NDArray[float32]) -> MotifBenchMeasurements: ...
def measure_ame(generated_catalytic_backbone: NDArray[float32], predicted_catalytic_backbone: NDArray[float32], generated_catalytic_heavy: NDArray[float32], predicted_catalytic_heavy: NDArray[float32], predicted_ligand: NDArray[float32], predicted_backbone: NDArray[float32]) -> AmeMeasurements: ...
def motifbench_1_0() -> CompatibilityProfile: ...
def ame_heavy_atom_1_0() -> CompatibilityProfile: ...
def read_evaluation_specification(path: PathLike[str]) -> EvaluationSpecification: ...

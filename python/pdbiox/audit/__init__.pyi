from collections.abc import Callable, Iterable
from typing import Any, final
from ..chem import RadiusSet
from ..query import AltlocPolicy, AnalysisPolicy, AssemblyChoice, MissingPolicy, ModelChoice, Namespace

@final
class HydrogenPolicy:
    ExplicitOnly: HydrogenPolicy
    IncludeInferred: HydrogenPolicy
    Exclude: HydrogenPolicy

@final
class EquivalencePolicy:
    @staticmethod
    def none() -> EquivalencePolicy: ...
    Ccd: EquivalencePolicy
    Explicit: EquivalencePolicy

@final
class SymmetryPolicy:
    @staticmethod
    def none() -> SymmetryPolicy: ...
    Crystallographic: SymmetryPolicy
    BiologicalAssembly: SymmetryPolicy

@final
class Precision:
    F32: Precision
    F64: Precision

@final
class PeriodicPolicy:
    @staticmethod
    def none() -> PeriodicPolicy: ...
    Pbc: PeriodicPolicy
    MinimumImage: PeriodicPolicy

@final
class AlignmentPolicy:
    @staticmethod
    def none() -> AlignmentPolicy: ...
    @staticmethod
    def explicit(selection: str) -> AlignmentPolicy: ...
    @staticmethod
    def global_( ) -> AlignmentPolicy: ...
    @staticmethod
    def local() -> AlignmentPolicy: ...

@final
class ContactDefinition:
    @staticmethod
    def distance_cutoff(tolerance: float) -> ContactDefinition: ...
    @staticmethod
    def surface_based(probe: float) -> ContactDefinition: ...

@final
class Tolerance:
    def __init__(self, relative: float, absolute: float) -> None: ...
    @staticmethod
    def standard() -> Tolerance: ...
    relative: float
    absolute: float

@final
class PolicyField:
    Assembly: PolicyField
    Model: PolicyField
    Altloc: PolicyField
    Identifiers: PolicyField
    MissingAtoms: PolicyField
    Hydrogens: PolicyField
    AtomEquivalence: PolicyField
    Symmetry: PolicyField
    Alignment: PolicyField
    Precision: PolicyField
    Periodic: PolicyField
    VdwRadii: PolicyField
    ContactDef: PolicyField
    FloatTolerance: PolicyField

@final
class PolicyValue:
    field: PolicyField
    value: str
    @staticmethod
    def named(field: PolicyField, value: str) -> PolicyValue: ...

@final
class PolicyDimension:
    @staticmethod
    def assembly(values: list[AssemblyChoice]) -> PolicyDimension: ...
    @staticmethod
    def model(values: list[ModelChoice]) -> PolicyDimension: ...
    @staticmethod
    def altloc(values: list[AltlocPolicy]) -> PolicyDimension: ...
    @staticmethod
    def identifiers(values: list[Namespace]) -> PolicyDimension: ...
    @staticmethod
    def missing_atoms(values: list[MissingPolicy]) -> PolicyDimension: ...
    @staticmethod
    def hydrogens(values: list[HydrogenPolicy]) -> PolicyDimension: ...
    @staticmethod
    def atom_equivalence(values: list[EquivalencePolicy]) -> PolicyDimension: ...
    @staticmethod
    def symmetry(values: list[SymmetryPolicy]) -> PolicyDimension: ...
    @staticmethod
    def alignment(values: list[AlignmentPolicy]) -> PolicyDimension: ...
    @staticmethod
    def precision(values: list[Precision]) -> PolicyDimension: ...
    @staticmethod
    def periodic(values: list[PeriodicPolicy]) -> PolicyDimension: ...
    @staticmethod
    def vdw_radii(values: list[RadiusSet]) -> PolicyDimension: ...
    @staticmethod
    def contact_def(values: list[ContactDefinition]) -> PolicyDimension: ...
    @staticmethod
    def float_tolerance(values: list[Tolerance]) -> PolicyDimension: ...
    field: PolicyField
    is_empty: bool
    def __len__(self) -> int: ...

@final
class PlanError:
    message: str
    field: PolicyField | None
    cost: int | None
    limit: int | None

@final
class PolicySpace:
    def __init__(self, policy: AnalysisPolicy | None = ...) -> None: ...
    def vary(self, dimension: PolicyDimension) -> PolicySpace: ...
    def with_max_runs(self, max_runs: int) -> PolicySpace: ...
    def cost(self) -> int: ...
    def plan(self) -> AuditPlan: ...

@final
class AuditPlan:
    cost: int
    fields: list[PolicyField]
    policies: list[AnalysisPolicy]

@final
class AuditRun:
    policy: AnalysisPolicy
    @property
    def result(self) -> Any: ...

@final
class SensitiveItem:
    item: str
    present_in: list[int]

@final
class DimensionSensitivity:
    field: PolicyField
    sensitive_items: list[str]
    mean_change: float

@final
class AuditReport:
    runs: list[AuditRun]
    stability: float
    sensitive_items: list[SensitiveItem]
    dimensions: list[DimensionSensitivity]

@final
class BatchDimension:
    field: PolicyField
    mean_change: float

@final
class BatchAudit:
    subjects: int
    mean_stability: float
    min_stability: float
    max_stability: float
    dimensions: list[BatchDimension]

def audit(plan: AuditPlan, analyse: Callable[[AnalysisPolicy], Any], items: Callable[[Any], Iterable[str]]) -> AuditReport: ...
def audit_batch(plan: AuditPlan, subjects: Iterable[Any], analyse: Callable[[Any, AnalysisPolicy], Any], items: Callable[[Any], Iterable[str]]) -> BatchAudit: ...

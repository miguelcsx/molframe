from typing import Generic, TypeVar, final

T = TypeVar("T")

@final
class Status:
    Complete: Status
    Partial: Status
    Ambiguous: Status
    Indeterminate: Status

@final
class Coverage:
    intended: int
    used: int
    missing: int
    ambiguous: int
    fraction: float | None
    is_complete: bool

@final
class ImpactEstimate:
    NoImpact: ImpactEstimate
    Low: ImpactEstimate
    Moderate: ImpactEstimate
    High: ImpactEstimate
    Unknown: ImpactEstimate

@final
class AssumptionSource:
    Explicit: AssumptionSource
    ProfileDefault: AssumptionSource
    Inferred: AssumptionSource

@final
class Assumption:
    field: str
    value: str
    source: AssumptionSource
    impact: ImpactEstimate
    is_silent_and_material: bool

@final
class Diagnostic:
    code: str
    message: str
    remedy: str
    span: tuple[int, int, int, int] | None

@final
class Provenance:
    pdbiox_version: str
    input_source: str
    input_fingerprint: str | None
    policy_fingerprint: str
    profile: str | None
    algorithm_name: str | None
    algorithm_version: str | None
    parameters: dict[str, bool | int | float | str]
    schema_version: str | None
    component_version: str | None
    timestamp: str | None
    fingerprint: str

@final
class Analysis(Generic[T]):
    value: T
    status: Status
    coverage: Coverage
    warnings: list[Diagnostic]
    assumptions: list[Assumption]
    provenance: Provenance
    is_usable: bool
    silent_assumptions: list[Assumption]

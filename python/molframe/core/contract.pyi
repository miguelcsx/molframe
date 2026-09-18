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
    def __init__(self, code: object) -> None: ...
    code: str
    message: str
    remedy: str
    span: tuple[int, int, int, int] | None
    code_value: object
    severity_value: object
    category_value: str | None
    field_value: str | None
    row: int | None
    context: list[object]
    def with_message(self, message: str) -> Diagnostic: ...
    def with_severity(self, severity: object) -> Diagnostic: ...
    def at(self, span: object) -> Diagnostic: ...
    def in_category(self, category: str) -> Diagnostic: ...
    def about_field(self, field: str) -> Diagnostic: ...
    def at_row(self, row: int) -> Diagnostic: ...
    def is_error(self, strictness: object) -> bool: ...
    def render(self, source: bytes | None = ..., *, origin: str | None = ..., color: bool = ...) -> object: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

@final
class Provenance:
    molframe_version: str
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
class ReexecutionEnvironment:
    def __init__(self, *, molframe_version: str | None = ..., schema_version: str | None = ..., component_version: str | None = ...) -> None: ...
    @staticmethod
    def current() -> ReexecutionEnvironment: ...
    molframe_version: str
    schema_version: str | None
    component_version: str | None

@final
class Reexecution:
    value: object
    provenance: Provenance

class ReexecutionError(Exception): ...

def reexecute_from_provenance(
    provenance: Provenance,
    input: bytes,
    environment: ReexecutionEnvironment,
    run: object,
) -> Reexecution: ...

@final
class AlgorithmId:
    def __init__(self, name: str, version: str) -> None: ...
    name: str
    version: str

@final
class DictionaryVersion:
    def __init__(self, version: str) -> None: ...
    value: str

@final
class ProfileId:
    @staticmethod
    def default() -> ProfileId: ...
    value: str
    def __str__(self) -> str: ...

@final
class Fingerprint:
    def __init__(self, data: bytes) -> None: ...
    @staticmethod
    def of(data: bytes) -> Fingerprint: ...
    value: int

@final
class ParameterValue:
    @staticmethod
    def boolean(value: bool) -> ParameterValue: ...
    @staticmethod
    def integer(value: int) -> ParameterValue: ...
    @staticmethod
    def float(value: float) -> ParameterValue: ...
    @staticmethod
    def text(value: str) -> ParameterValue: ...
    kind: str
    bool_value: bool | None
    integer_value: int | None
    float_value: float | None
    text_value: str | None
    def to_python(self) -> bool | int | float | str: ...

@final
class AnalysisParameters:
    def __init__(self, values: dict[str, ParameterValue] | None = ...) -> None: ...
    def __len__(self) -> int: ...
    def __contains__(self, name: str) -> bool: ...
    def get(self, name: str) -> ParameterValue | None: ...
    def __getitem__(self, name: str) -> ParameterValue: ...
    def set(self, name: str, value: ParameterValue) -> None: ...
    def items(self) -> list[tuple[str, ParameterValue]]: ...
    def to_dict(self) -> dict[str, bool | int | float | str]: ...

@final
class SourceRef:
    @staticmethod
    def none() -> SourceRef: ...
    @staticmethod
    def path(path: str) -> SourceRef: ...
    @staticmethod
    def url(url: str) -> SourceRef: ...
    @staticmethod
    def memory() -> SourceRef: ...
    kind: str
    path_value: str | None
    def __str__(self) -> str: ...

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

from typing import final
from . import CifValue, Document

@final
class ModelRow:
    values: list[CifValue]
    def value(self, index: int) -> CifValue | None: ...

@final
class ModelCategory:
    name: str
    items: list[str]
    rows: list[ModelRow]
    def value(self, item: str, row: int) -> CifValue | None: ...

@final
class MetricDefinition:
    id: str
    name: str | None
    metric_type: str
    mode: str
    software_group_id: str | None

@final
class GlobalMetric:
    model_id: str
    metric_id: str
    value: float

@final
class LocalMetric:
    model_id: str
    chain_id: str
    sequence_id: int
    component_id: str | None
    metric_id: str
    value: float

@final
class PairwiseMetric:
    model_id: str
    first_chain_id: str
    first_sequence_id: int
    second_chain_id: str
    second_sequence_id: int
    metric_id: str
    value: float

@final
class QualityMetrics:
    definitions: list[MetricDefinition]
    global_metrics: list[GlobalMetric]
    local: list[LocalMetric]
    pairwise: list[PairwiseMetric]
    def plddt(self) -> list[LocalMetric]: ...
    def pae(self) -> list[PairwiseMetric]: ...
    def ptm(self) -> list[GlobalMetric]: ...

@final
class ModelDescription:
    id: str
    assembly_id: str | None
    name: str | None
    model_type: str | None

@final
class Target:
    entity_id: str
    data_id: str | None
    origin: str | None

@final
class Template:
    id: str
    target_chain_id: str | None
    name: str | None
    origin: str | None
    entity_type: str | None

@final
class ProtocolStep:
    protocol_id: str
    step_id: str
    method_type: str
    name: str | None
    details: str | None
    software_group_id: str | None

@final
class SoftwareGroup:
    group_id: str
    software_id: str
    parameter_group_id: str | None

@final
class ModelCif:
    is_empty: bool
    categories: list[ModelCategory]
    models: list[ModelDescription]
    targets: list[Target]
    templates: list[Template]
    protocol_steps: list[ProtocolStep]
    software_groups: list[SoftwareGroup]
    quality: QualityMetrics
    def confidence(self) -> QualityMetrics: ...

def lower_model_cif(document: Document) -> ModelCif: ...
def lower(document: Document) -> ModelCif: ...
def modelcif_write_canonical(structure: object, model_cif: ModelCif, options: object | None = ...) -> str: ...
def modelcif_write_canonical_with_options(structure: object, model_cif: ModelCif, options: object) -> str: ...
def write_canonical(structure: object, model_cif: ModelCif, options: object | None = ...) -> str: ...
def write_canonical_with_options(structure: object, model_cif: ModelCif, options: object) -> str: ...
def modelcif_attach(structure: object, model_cif: ModelCif) -> object: ...
MODEL_CIF_EXTENSION: str

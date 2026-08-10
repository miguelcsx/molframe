from typing import final
from . import Structure
from .contract import Analysis
from .query import AnalysisPolicy, Selection

@final
class BFactorOutlier:
    atom: int; value: float; z_score: float
@final
class BFactorDistribution:
    intended: int; assessed: int; mean: float; variance: float
    standard_deviation: float; minimum: float; median: float; maximum: float
    outliers: list[BFactorOutlier]
@final
class TlsModel:
    def __init__(self, origin: list[float], translation: list[list[float]], libration: list[list[float]], screw: list[list[float]]) -> None: ...
@final
class TlsGroup:
    def __init__(self, id: str, atoms: Selection, model: TlsModel) -> None: ...
@final
class TlsBFactorFlag:
    group: str; atom: int; observed: float; predicted: float; deviation: float
@final
class TlsBFactorReport:
    intended: int; assessed: int; flags: list[TlsBFactorFlag]

def b_factor_distribution(structure: Structure, selection: Selection, outlier_standard_deviations: float) -> BFactorDistribution: ...
def analyse_b_factor_distribution(structure: Structure, selection: Selection, outlier_standard_deviations: float, policy: AnalysisPolicy) -> Analysis[BFactorDistribution]: ...
def tls_b_factor_consistency(structure: Structure, groups: list[TlsGroup], maximum_absolute_deviation: float, symmetry_tolerance: float) -> TlsBFactorReport: ...
def analyse_tls_b_factor_consistency(structure: Structure, groups: list[TlsGroup], maximum_absolute_deviation: float, symmetry_tolerance: float, policy: AnalysisPolicy) -> Analysis[TlsBFactorReport]: ...

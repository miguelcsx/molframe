from typing import final
from numpy import float64
from numpy.typing import NDArray

@final
class XvgSeries:
    title: str | None
    abscissa_label: str | None
    ordinate_label: str | None
    legends: list[str]
    abscissa: NDArray[float64]
    def __len__(self) -> int: ...
    def is_empty(self) -> bool: ...
    def column(self, index: int) -> NDArray[float64]: ...
    def column_by_legend(self, legend: str) -> NDArray[float64] | None: ...
    def sample(self, column: int, abscissa: float) -> float | None: ...

def read_xvg(data: bytes | str) -> XvgSeries: ...

__all__: list[str]

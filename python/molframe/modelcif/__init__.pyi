from os import PathLike
from typing import Literal

import numpy as np

from .. import Structure

class PaeMatrix:
    """Square residue-pair aligned error, stored in its narrowest exact layout."""

    @property
    def size(self) -> int:
        """Residues per side."""
    @property
    def dtype(self) -> Literal["uint8", "uint16", "float16", "float32"]:
        """NumPy dtype name of the stored cells."""
    @property
    def denominator(self) -> int:
        """Divide fixed-point cells by this to get values; 1 for floating point."""
    @property
    def nbytes(self) -> int:
        """Bytes the cell buffer occupies."""
    @property
    def is_exact(self) -> bool:
        """Whether every source value is reproduced exactly."""
    @property
    def max_quantisation_error(self) -> float:
        """Largest absolute difference between a source value and its stored value."""
    @property
    def declared_maximum(self) -> float | None:
        """Maximum error the producing model declared, when it stated one."""
    def array(self) -> np.ndarray:
        """Read-only ``(size, size)`` zero-copy view of the stored cells."""
    def get(self, row: int, column: int) -> float:
        """Return the value at ``(row, column)`` as the source stated it."""
    def block_mean(self, rows: tuple[int, int], columns: tuple[int, int]) -> float | None:
        """Mean of a rectangle given as half-open ranges."""
    def mean(self) -> float | None:
        """Mean over every cell."""
    def row(self, row: int) -> list[float]:
        """One row decoded to float32 (a copy)."""

def pae(structure: Structure) -> PaeMatrix | None:
    """Return the dense PAE matrix of a structure read from a ModelCIF file, if any."""

def read_pae_json(
    path: str | PathLike[str],
    *,
    quantise: Literal["exact", "uint8", "uint16", "float16", "float32"] | None = None,
    max_bytes: int | None = None,
) -> PaeMatrix:
    """Import an AlphaFold-style PAE JSON file as a dense matrix."""

PLDDT_UNSCORED: int

def plddt_bands(structure: Structure) -> np.ndarray | None:
    """Return each atom's confidence band as uint8, or None without pLDDT.

    Codes index ``plddt_band_names()``; ``PLDDT_UNSCORED`` marks atoms with no score.
    """

def plddt_band_names() -> list[str]:
    """Return the band names, least to most confident."""

def plddt_band_colors() -> np.ndarray:
    """Return the conventional sRGB colour of each band as a ``(4, 3)`` uint8 array."""

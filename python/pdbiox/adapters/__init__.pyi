from typing import final
from ..core import Structure
from ..query import Namespace
from .._io_types import BondOrder

class DownloadError(RuntimeError): ...
class TopologyExportError(ValueError): ...

@final
class DownloadOptions:
    def __init__(self, *, max_bytes: int, timeout_seconds: float, redirect_limit: int = ...) -> None: ...
    max_bytes: int
    timeout_seconds: float
    redirect_limit: int

@final
class VerifiedDownload:
    bytes: bytes
    sha256: str

@final
class ExportChain:
    id: str
    residues: tuple[int, int]

@final
class ExportResidue:
    name: str
    number: int | None
    insertion_code: str | None
    is_heterogen: bool
    chain: int
    atoms: tuple[int, int]

@final
class ExportAtom:
    name: str
    atomic_number: int
    element_symbol: str
    mass: float
    serial: int | None
    formal_charge: int | None
    occupancy: float | None
    b_factor: float | None
    alternate_location: str | None
    residue: int
    position: list[float]

@final
class ExportBond:
    atom_a: int
    atom_b: int
    order: BondOrder

@final
class TopologyExport:
    @classmethod
    def from_model(cls, structure: Structure, model: int, namespace: Namespace) -> TopologyExport: ...
    chains: list[ExportChain]
    residues: list[ExportResidue]
    atoms: list[ExportAtom]
    bonds: list[ExportBond]

def fetch_verified(url: str, expected_sha256: str, options: DownloadOptions) -> VerifiedDownload: ...

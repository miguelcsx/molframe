from typing import final
from ..core import Structure
from ..query import Namespace
from .._io_types import BondOrder

MISSING_STRING: int

class DownloadError(RuntimeError): ...
class TopologyBatchError(ValueError): ...

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
class TopologyBatch:
    @classmethod
    def from_model(cls, structure: Structure, model: int, namespace: Namespace) -> TopologyBatch: ...
    strings: list[str]
    chain_ids: list[int]
    chain_residue_offsets: list[int]
    residue_names: list[int]
    residue_numbers: list[int]
    residue_number_validity: list[bool]
    residue_insertion_codes: list[int]
    residue_is_heterogen: list[bool]
    residue_chain: list[int]
    residue_atom_offsets: list[int]
    atom_names: list[int]
    atomic_numbers: list[int]
    masses: list[float]
    atom_serials: list[int]
    atom_serial_validity: list[bool]
    formal_charges: list[int]
    formal_charge_validity: list[bool]
    occupancies: list[float]
    occupancy_validity: list[bool]
    b_factors: list[float]
    b_factor_validity: list[bool]
    atom_alternate_locations: list[int]
    atom_residue: list[int]
    position_x: list[float]
    position_y: list[float]
    position_z: list[float]
    bond_atom_a: list[int]
    bond_atom_b: list[int]
    bond_orders: list[BondOrder]
    def transfer_to_structure(self) -> Structure: ...

def fetch_verified(url: str, expected_sha256: str, options: DownloadOptions) -> VerifiedDownload: ...

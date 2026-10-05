from collections.abc import Mapping
from os import PathLike
from typing import Literal

from .. import Structure

def to_mmcif(structure: Structure) -> str: ...
def to_pdb(
    structure: Structure,
    *,
    hybrid36: bool = False,
    chain_map: Mapping[str, str] | None = None,
) -> str: ...
def to_bcif(
    structure: Structure,
    *,
    block_id: str | None = None,
    generated_connections: bool = False,
    transport: bool = False,
) -> bytes:
    """Encode BinaryCIF with explicit identities; transport regenerates label namespaces."""

def write(
    structure: Structure,
    path: str | PathLike[str],
    *,
    format: Literal["mmcif", "bcif", "pdb", "pqr", "pdbqt", "mmtf", "sdf", "mol2"] | None = None,
    hybrid36: bool = False,
    chain_map: Mapping[str, str] | None = None,
    memory_limit: int | None = None,
) -> None: ...

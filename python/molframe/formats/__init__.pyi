from collections.abc import Mapping
from os import PathLike
from typing import Literal

from .. import Structure

def to_mmcif(
    structure: Structure,
    *,
    block_id: str | None = None,
    generated_connections: bool = True,
) -> str:
    """Write mmCIF keeping space group, resolution, assemblies, NCS, secondary structure."""

def preserve_mmcif(source: str | PathLike[str], destination: str | PathLike[str]) -> None:
    """Copy an mmCIF document category for category, losing nothing."""

def to_pdb(
    structure: Structure,
    *,
    hybrid36: bool = False,
    chain_map: Mapping[str, str] | None = None,
    namespace: Literal["auto", "auth", "label"] = "auto",
) -> str:
    """Write PDB; ``auto`` uses author identifiers when every residue has them, else label."""

def to_bcif(
    structure: Structure,
    *,
    block_id: str | None = None,
    generated_connections: bool = True,
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
    namespace: Literal["auto", "auth", "label"] = "auto",
    block_id: str | None = None,
    generated_connections: bool = True,
    memory_limit: int | None = None,
) -> None: ...

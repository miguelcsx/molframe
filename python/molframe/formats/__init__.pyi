from __future__ import annotations
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
def to_bcif(structure: Structure) -> bytes: ...
def write(
    structure: Structure,
    path: str | PathLike[str],
    *,
    format: Literal["mmcif", "cif", "pdb", "ent", "bcif"] | None = None,
) -> None: ...

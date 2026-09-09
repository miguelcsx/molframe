from typing import Any
from ..core import Structure
from .._native import HydrogenBond, HydrogenBondOptions

def hydrogen_bonds(structure: Structure, options: Any, context: object | None = None) -> list[Any]: ...

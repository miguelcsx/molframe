from typing import final

@final
class AmberSection:
    format: str
    values: list[str]

@final
class AmberTopologyAtom:
    name: str
    charge: float
    mass: float
    type_index: int
    atom_type: str
    atomic_number: int | None
    residue: int

@final
class AmberTopologyResidue:
    name: str
    atom_start: int
    atom_end: int

@final
class AmberTopologyBond:
    atoms: tuple[int, int]
    parameter: int
    includes_hydrogen: bool

@final
class AmberTopologyAngle:
    atoms: tuple[int, int, int]
    parameter: int
    includes_hydrogen: bool

@final
class AmberTopologyDihedral:
    atoms: tuple[int, int, int, int]
    parameter: int
    improper: bool
    ignore_end_group: bool
    includes_hydrogen: bool

@final
class AmberTopology:
    version: str | None
    atoms: list[AmberTopologyAtom]
    residues: list[AmberTopologyResidue]
    bonds: list[AmberTopologyBond]
    angles: list[AmberTopologyAngle]
    dihedrals: list[AmberTopologyDihedral]
    sections: dict[str, AmberSection]

class AmberTopologyError(Exception): ...
def parse_amber_topology(text: str) -> AmberTopology: ...

__all__: list[str]

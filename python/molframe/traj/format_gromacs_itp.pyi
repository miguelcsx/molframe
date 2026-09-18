from typing import final

@final
class GromacsMoleculeType:
    name: str
    exclusions: int

@final
class GromacsItpAtom:
    number: int
    atom_type: str
    residue_number: int
    residue_name: str
    atom_name: str
    charge_group: int
    charge: float
    mass: float | None
    state_b: list[str]

@final
class GromacsInteraction:
    atoms: list[int]
    function: int
    parameters: list[str]

@final
class GromacsItp:
    directives: list[str]
    molecule_type: GromacsMoleculeType | None
    atoms: list[GromacsItpAtom]
    bonds: list[GromacsInteraction]
    pairs: list[GromacsInteraction]
    angles: list[GromacsInteraction]
    dihedrals: list[GromacsInteraction]
    constraints: list[GromacsInteraction]
    exclusions: list[list[int]]
    other_sections: dict[str, list[str]]

class GromacsItpError(Exception): ...
def parse_gromacs_itp(text: str) -> GromacsItp: ...

__all__: list[str]

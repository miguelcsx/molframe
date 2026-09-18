from typing import final

@final
class PsfAtom:
    segment: str
    residue_id: str
    residue_name: str
    atom_name: str
    atom_type: str
    charge: float
    mass: float

@final
class PsfTopology:
    titles: list[str]
    atoms: list[PsfAtom]
    bonds: list[tuple[int, int]]
    angles: list[tuple[int, int, int]]
    dihedrals: list[tuple[int, int, int, int]]
    impropers: list[tuple[int, int, int, int]]
    donors: list[tuple[int, int]]
    acceptors: list[tuple[int, int]]
    cross_terms: list[tuple[int, int, int, int, int, int, int, int]]

class PsfError(Exception): ...
def parse_psf(text: str) -> PsfTopology: ...
def write_psf(topology: PsfTopology) -> str: ...

__all__: list[str]

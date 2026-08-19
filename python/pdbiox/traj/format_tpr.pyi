class TprError(Exception): ...
class TprHeader:
    producer: str
    format_version: int
    generation: int
    precision: int
    atom_count: int
    tag: str
    has_coordinates: bool
class TprResidue:
    name: str
    number: int
    molecule: int
    molecule_type: str
class TprAtom:
    index: int
    name: str
    atom_type: str
    residue: int
    mass: float
    charge: float
    atomic_number: int | None
class TprBond:
    atom_a: int
    atom_b: int
class TprTopology:
    header: TprHeader
    residues: list[TprResidue]
    atoms: list[TprAtom]
    bonds: list[TprBond]
def parse_tpr(bytes: bytes) -> TprTopology: ...

__all__: list[str]

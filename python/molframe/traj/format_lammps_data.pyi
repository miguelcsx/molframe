class LammpsDataError(Exception): ...
class LammpsAtomStyle:
    Atomic: LammpsAtomStyle
    Charge: LammpsAtomStyle
    Molecular: LammpsAtomStyle
    Full: LammpsAtomStyle
class LammpsDataAtom:
    id: int
    molecule: int | None
    atom_type: int
    charge: float | None
    position: tuple[float, float, float]
    image: tuple[int, int, int] | None
class LammpsInteraction:
    id: int
    interaction_type: int
    atoms: list[int]
class LammpsDataCell:
    lower: tuple[float, float, float]
    upper: tuple[float, float, float]
    tilt: tuple[float, float, float]
class LammpsData:
    title: str
    atom_style: LammpsAtomStyle
    cell: LammpsDataCell | None
    masses: dict[int, float]
    atoms: list[LammpsDataAtom]
    bonds: list[LammpsInteraction]
    angles: list[LammpsInteraction]
    dihedrals: list[LammpsInteraction]
    impropers: list[LammpsInteraction]
    other_sections: dict[str, list[str]]
def parse_lammps_data(source: str) -> LammpsData: ...

__all__: list[str]

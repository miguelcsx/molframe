from typing import final
from os import PathLike

@final
class ComponentDictionary:
    def __init__(self, path: str | PathLike[str], version: str) -> None: ...
    version: str
    def peoe_charges(self, component_id: str, options: PeoeOptions) -> list[float]: ...

@final
class RadiusSet:
    Bondi: RadiusSet
    AmberUnited: RadiusSet
    Charmm: RadiusSet
    Alvarez: RadiusSet
    name: str
    version: str

@final
class ElementProperties:
    atomic_weight: float
    covalent_radius: float | None
    electronegativity: float | None
    valence_electrons: int
    period: int
    group: int | None

@final
class IonicSpin:
    Unspecified: IonicSpin
    High: IonicSpin
    Low: IonicSpin

@final
class IonicRadius:
    charge: int
    coordination: str
    spin: IonicSpin
    ionic_radius: float
    crystal_radius: float
    most_reliable: bool | None

@final
class Element:
    def __init__(self, symbol: str) -> None: ...
    @staticmethod
    def from_atomic_number(atomic_number: int) -> Element: ...
    symbol: str
    atomic_number: int
    properties: ElementProperties | None
    def vdw_radius(self, set: RadiusSet) -> float | None: ...
    def ionic_radii(self) -> list[IonicRadius]: ...
    def __repr__(self) -> str: ...

@final
class PeoeOptions:
    def __init__(self, passes: int, initial_damping: float, damping_factor: float, minimum_electronegativity_difference: float, profile: PeoeParameterProfile) -> None: ...
    @staticmethod
    def standard() -> PeoeOptions: ...

@final
class PeoeParameterProfile:
    GasteigerMarsili: PeoeParameterProfile

@final
class PeoeAtomType:
    H: PeoeAtomType
    CSp3: PeoeAtomType
    CSp2: PeoeAtomType
    CSp: PeoeAtomType
    NSp3: PeoeAtomType
    NSp2: PeoeAtomType
    NSp: PeoeAtomType
    OSp3: PeoeAtomType
    OSp2: PeoeAtomType
    FSp3: PeoeAtomType
    ClSp3: PeoeAtomType
    BrSp3: PeoeAtomType
    ISp3: PeoeAtomType
    SSp3: PeoeAtomType
    SO: PeoeAtomType
    SO2: PeoeAtomType
    SSp2: PeoeAtomType
    PSp3: PeoeAtomType
    PSp2: PeoeAtomType
    SiSp3: PeoeAtomType
    SiSp2: PeoeAtomType
    SiSp: PeoeAtomType
    BSp3: PeoeAtomType
    BSp2: PeoeAtomType
    BeSp3: PeoeAtomType
    BeSp2: PeoeAtomType
    MgSp3: PeoeAtomType
    MgSp2: PeoeAtomType
    MgSp: PeoeAtomType
    AlSp3: PeoeAtomType
    AlSp2: PeoeAtomType

@final
class PeoeAtom:
    def __init__(self, atom_type: PeoeAtomType, formal_charge: float) -> None: ...

@final
class PeoeBond:
    def __init__(self, atom_a: int, atom_b: int) -> None: ...

def peoe_charges(atoms: list[PeoeAtom], bonds: list[PeoeBond], options: PeoeOptions) -> list[float]: ...

from typing import final
from os import PathLike
from . import BondInference, BondInferenceReport, Structure

DEFAULT_BOND_RADIUS_SCALE: float
DEFAULT_MINIMUM_BOND_DISTANCE: float

@final
class InferBonds:
    def __init__(self, *, options: BondInference | None = ...) -> None: ...
    options: BondInference
    def execute(self, structure: Structure) -> BondInferenceReport: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> InferBonds: ...
    def __repr__(self) -> str: ...

@final
class ComponentDictionary:
    def __init__(self, path: str | PathLike[str], version: str) -> None: ...
    version: str
    def get(self, component_id: str) -> Component | None: ...
    def equivalence_classes(self, component_id: str) -> EquivalenceClasses | None: ...
    def coverage(self, structure: object, policy: object) -> ComponentCoverage: ...
    def coverage_with_policy(self, structure: object, policy: object) -> ComponentCoverage: ...
    def side_chain_definition(self, component_id: str, nitrogen: str, alpha_carbon: str, side_chain_atoms: list[str]) -> list[str] | None: ...
    def peoe_charges(self, component_id: str, options: PeoeOptions) -> list[float]: ...
    def apply_chemistry(self, structure: object, policy: PolymerLinkPolicy | None = None) -> ChemistryReport: ...
    def apply_polymer_role_profile(self, structure: object, profile: PolymerRoleProfile) -> PolymerRoleReport: ...

CifProvider = ComponentDictionary
ComponentProvider = ComponentDictionary

@final
class MemoryProvider:
    def __init__(self, version: str, components: list[Component]) -> None: ...
    version: str
    def get(self, component_id: str) -> Component | None: ...

@final
class RadiusTable:
    name: str
    version: str

class PeoeError(ValueError): ...
class AutomorphismLimit(ValueError): ...

@final
class ComponentKind:
    AminoAcid: ComponentKind
    Nucleotide: ComponentKind
    Saccharide: ComponentKind
    Lipid: ComponentKind
    NonPolymer: ComponentKind
    Solvent: ComponentKind
    Ion: ComponentKind
    Unknown: ComponentKind

@final
class StereoConfiguration:
    R: StereoConfiguration
    S: StereoConfiguration
    Mixed: StereoConfiguration

@final
class ComponentAtom:
    name: str
    alternate_name: str | None
    element: Element
    charge: int
    aromatic: bool
    leaving: bool
    stereo: StereoConfiguration | None

@final
class ComponentBond:
    atom_a: str
    atom_b: str
    order: object
    aromatic: bool
    stereo: StereoConfiguration | None

@final
class Component:
    id: str
    name: str
    kind: ComponentKind
    parent: str | None
    one_letter_code: str | None
    formula: str | None
    atoms: list[ComponentAtom]
    bonds: list[ComponentBond]
    ideal_coordinates: list[tuple[float, float, float]] | None
    model_coordinates: list[tuple[float, float, float]] | None
    def atom(self, name: str) -> ComponentAtom | None: ...

@final
class ComponentCoverage:
    coverage: object
    findings: list[str]
    dictionary_version: str

@final
class EquivalenceClasses:
    classes: list[list[int]]
    def class_of(self, atom: int) -> int | None: ...
    def equivalent(self, atom_a: int, atom_b: int) -> bool: ...

@final
class EquivalenceCache:
    def __init__(self, version: str) -> None: ...
    version: str
    len: int
    is_empty: bool
    def get(self, component: Component) -> EquivalenceClasses: ...

@final
class ChemistryReport:
    structure: object
    findings: list[str]
    dictionary_version: str
    polymer_link_policy: PolymerLinkPolicy

@final
class PolymerLinkRule:
    def __init__(self, left_kind: ComponentKind, left_atom: str, right_kind: ComponentKind, right_atom: str) -> None: ...
    left_kind: ComponentKind
    left_atom: str
    right_kind: ComponentKind
    right_atom: str

@final
class PolymerLinkPolicy:
    @staticmethod
    def disabled() -> PolymerLinkPolicy: ...
    @staticmethod
    def explicit(angstrom: float, rules: list[PolymerLinkRule]) -> PolymerLinkPolicy: ...
    is_disabled: bool

@final
class PolymerRoleRule:
    def __init__(self, component_id: str | None, component_kind: ComponentKind | None, atom_name: str, role: PolymerAtomRole) -> None: ...
    component_id: str | None
    component_kind: ComponentKind | None
    atom_name: str
    role: PolymerAtomRole

@final
class PolymerRoleProfile:
    def __init__(self, id: str, rules: list[PolymerRoleRule]) -> None: ...
    id: str
    rules: list[PolymerRoleRule]

@final
class PolymerRoleReport:
    structure: object
    unresolved_components: list[str]
    dictionary_version: str
    profile_id: str

@final
class PolymerAtomRole:
    @staticmethod
    def from_code(code: int) -> PolymerAtomRole | None: ...
    code: int
    def intersects(self, required: PolymerAtomRole) -> bool: ...
    def union(self, other: PolymerAtomRole) -> PolymerAtomRole: ...
    unknown: PolymerAtomRole
    protein_backbone: PolymerAtomRole
    nucleic_backbone: PolymerAtomRole

@final
class RadiusSet:
    Bondi: RadiusSet
    AmberUnited: RadiusSet
    Charmm: RadiusSet
    Alvarez: RadiusSet
    name: str
    version: str
    def table(self) -> RadiusTable: ...

RadiiSet = RadiusSet

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
    def from_symbol(symbol: str) -> Element | None: ...
    @staticmethod
    def from_atomic_number(atomic_number: int) -> Element: ...
    @staticmethod
    def infer_from_name(name: str) -> Element: ...
    @staticmethod
    def infer_from_pdb_atom_name(field: str) -> Element: ...
    UNKNOWN: Element
    HYDROGEN: Element
    CARBON: Element
    NITROGEN: Element
    OXYGEN: Element
    FLUORINE: Element
    PHOSPHORUS: Element
    SULFUR: Element
    CHLORINE: Element
    CALCIUM: Element
    IRON: Element
    ZINC: Element
    SELENIUM: Element
    MAX_ATOMIC_NUMBER: int
    symbol: str
    atomic_number: int
    properties: ElementProperties | None
    def vdw_radius(self, set: RadiusSet) -> float | None: ...
    def ionic_radii(self) -> list[IonicRadius]: ...
    def is_unknown(self) -> bool: ...
    def is_hydrogen(self) -> bool: ...
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

def read_ccd(path: str | PathLike[str], version: str) -> tuple[ComponentDictionary, list[str]]: ...
def element_properties(element: Element) -> ElementProperties | None: ...
def vdw_radius(element: Element, set: RadiusSet) -> float | None: ...
def ionic_radii(element: Element) -> list[IonicRadius]: ...
def component_coverage(structure: object, provider: CifProvider | MemoryProvider, policy: object) -> ComponentCoverage: ...
def component_peoe_charges(component: Component, options: PeoeOptions) -> list[float]: ...
def apply_component_chemistry(structure: object, provider: CifProvider | MemoryProvider, policy: PolymerLinkPolicy | None = ...) -> ChemistryReport: ...
def apply_polymer_role_profile(structure: object, provider: CifProvider | MemoryProvider, profile: PolymerRoleProfile) -> PolymerRoleReport: ...
def equivalence_classes(component: Component) -> EquivalenceClasses: ...
def automorphisms(component: Component, limit: int) -> list[list[int]]: ...

@final
class SideChainRoles:
    def __init__(self, nitrogen: str, alpha_carbon: str, side_chain_atoms: list[str]) -> None: ...
    nitrogen: str
    alpha_carbon: str
    side_chain_atoms: list[str]

@final
class SideChainDefinition:
    atoms: list[str]
    def torsion_count(self) -> int: ...

def side_chain_definition(component: Component, roles: SideChainRoles) -> SideChainDefinition | None: ...

class SmartsError(ValueError): ...
class SmartsDataError(ValueError): ...

@final
class SmartsMatch:
    atom_indices: list[int]

@final
class SmartsPattern:
    def __init__(self, text: str) -> None: ...
    @staticmethod
    def parse(text: str) -> SmartsPattern: ...
    def find_matches(self, component: Component) -> list[SmartsMatch]: ...
    def find_structure_matches(self, structure: object) -> list[SmartsMatch]: ...
    def matches(self, component: Component) -> bool: ...
    def atom_count(self) -> int: ...

def parse_smarts(text: str) -> SmartsPattern: ...

class MolError(ValueError): ...

class Mol2Error(ValueError): ...

@final
class MolAtom:
    def __init__(self, element: Element, position: list[float]) -> None: ...
    element: Element
    position: tuple[float, float, float]

@final
class MolBond:
    def __init__(self, first: int, second: int, order: int) -> None: ...
    first: int
    second: int
    order: int

@final
class Molecule:
    def __init__(self, atoms: list[MolAtom], bonds: list[MolBond]) -> None: ...
    atoms: list[MolAtom]
    bonds: list[MolBond]

@final
class MolVersion:
    V2000: MolVersion
    V3000: MolVersion

@final
class MolAtomMetadata:
    def __init__(self, *, formal_charge: int | None = None, isotope: int | None = None, stereo_parity: int | None = None) -> None: ...
    formal_charge: int | None
    isotope: int | None
    stereo_parity: int | None

@final
class MolBondMetadata:
    def __init__(self, *, stereo: int | None = None) -> None: ...
    stereo: int | None

@final
class SdfProperty:
    def __init__(self, name: str, value: str) -> None: ...
    name: str
    value: str

@final
class MolRecord:
    def __init__(self, name: str, program: str, comment: str, version: MolVersion, molecule: Molecule, atom_metadata: list[MolAtomMetadata], bond_metadata: list[MolBondMetadata], properties: list[SdfProperty]) -> None: ...
    name: str
    program: str
    comment: str
    version: MolVersion
    molecule: Molecule
    atom_metadata: list[MolAtomMetadata]
    bond_metadata: list[MolBondMetadata]
    properties: list[SdfProperty]

@final
class Mol2AtomMetadata:
    def __init__(self, id: int, name: str, atom_type: str, *, substructure_id: int | None = None, substructure_name: str | None = None, charge: float | None = None, status_bits: str | None = None) -> None: ...
    id: int
    name: str
    atom_type: str
    substructure_id: int | None
    substructure_name: str | None
    charge: float | None
    status_bits: str | None

@final
class Mol2BondMetadata:
    def __init__(self, id: int, bond_type: str, *, status_bits: str | None = None) -> None: ...
    id: int
    bond_type: str
    status_bits: str | None

@final
class Mol2Section:
    def __init__(self, name: str, lines: list[str]) -> None: ...
    name: str
    lines: list[str]

@final
class Mol2Record:
    def __init__(self, name: str, molecule_type: str, charge_type: str, additional_counts: list[int], status_bits: str | None, comment: str | None, molecule: Molecule, atom_metadata: list[Mol2AtomMetadata], bond_metadata: list[Mol2BondMetadata], extra_sections: list[Mol2Section]) -> None: ...
    name: str
    molecule_type: str
    charge_type: str
    additional_counts: list[int]
    status_bits: str | None
    comment: str | None
    molecule: Molecule
    atom_metadata: list[Mol2AtomMetadata]
    bond_metadata: list[Mol2BondMetadata]
    extra_sections: list[Mol2Section]

def parse_mol_record(text: str) -> MolRecord: ...
def parse_sdf_records(text: str) -> list[MolRecord]: ...
def write_mol(record: MolRecord) -> str: ...
def write_sdf(records: list[MolRecord]) -> str: ...
def parse_mol2_record(text: str) -> Mol2Record: ...
def write_mol2(record: Mol2Record) -> str: ...

from collections.abc import Mapping, Sequence
from os import PathLike
from typing import Literal, Protocol, final

from .. import Structure

type _RadiusSet = Literal["bondi", "amber_united", "charmm", "alvarez"]

class Float32Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class UInt32Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

@final
class ElementProperties:
    @property
    def symbol(self) -> str: ...
    @property
    def atomic_number(self) -> int: ...
    @property
    def atomic_weight(self) -> float: ...
    @property
    def covalent_radius(self) -> float | None: ...
    @property
    def electronegativity(self) -> float | None: ...
    @property
    def valence_electrons(self) -> int: ...
    @property
    def period(self) -> int: ...
    @property
    def group(self) -> int | None: ...

def element(symbol: str) -> ElementProperties: ...
def vdw_radius(symbol: str, *, radii: _RadiusSet = "bondi") -> float | None: ...
def vdw_radii(structure: Structure, *, radii: _RadiusSet = "bondi") -> Float32Array: ...
def annotate(
    structure: Structure,
    components: str | PathLike[str],
    *,
    version: str = "unversioned",
) -> Structure: ...

@final
class PolymerRoleRule:
    def __init__(
        self,
        atom_name: str,
        role: int,
        *,
        component_id: str | None = None,
        component_kind: int | None = None,
    ) -> None: ...
    @property
    def atom_name(self) -> str: ...
    @property
    def role(self) -> int: ...
    @property
    def component_id(self) -> str | None: ...
    @property
    def component_kind(self) -> int | None: ...

@final
class PolymerRoleReport:
    @property
    def structure(self) -> Structure: ...
    @property
    def unresolved_components(self) -> list[str]: ...
    @property
    def dictionary_version(self) -> str: ...
    @property
    def profile_id(self) -> str: ...

def apply_polymer_role_profile(
    structure: Structure,
    components: str | PathLike[str],
    rules: list[PolymerRoleRule],
    *,
    profile_id: str,
    version: str = "unversioned",
) -> PolymerRoleReport: ...

@final
class PartialCharges:
    @property
    def values(self) -> list[float]: ...
    @property
    def source(self) -> Literal["file", "peoe"]: ...
    @property
    def dictionary_version(self) -> str | None: ...
    @property
    def parameter_profile(self) -> str | None: ...

def partial_charges(
    structure: Structure,
    components: str | PathLike[str] | None = None,
    *,
    version: str = "unversioned",
) -> PartialCharges: ...

type _SnfgShape = Literal[
    "filled_sphere",
    "filled_cube",
    "crossed_cube",
    "divided_diamond",
    "filled_cone",
    "divided_cone",
    "flat_box",
    "filled_star",
    "filled_diamond",
    "flat_diamond",
    "flat_hexagon",
    "pentagon",
    "diamond_prism",
    "pentagonal_prism",
    "hexagonal_prism",
    "heptagonal_prism",
]
type _BondProvenance = Literal["file", "chemical_component_dictionary", "inferred_distance", "user"]

@final
class SnfgSymbol:
    @property
    def abbreviation(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def color(self) -> int: ...
    @property
    def shape(self) -> _SnfgShape: ...
    @property
    def secondary_color(self) -> int | None: ...

@final
class RingGeometry:
    @property
    def center(self) -> tuple[float, float, float]: ...
    @property
    def normal(self) -> tuple[float, float, float]: ...
    @property
    def anomeric_direction(self) -> tuple[float, float, float] | None: ...

@final
class Monosaccharide:
    @property
    def residue(self) -> int: ...
    @property
    def ring_atoms(self) -> list[int]: ...
    @property
    def anomeric_atom(self) -> int | None: ...
    @property
    def symbol(self) -> SnfgSymbol: ...
    @property
    def geometry(self) -> RingGeometry | None: ...

@final
class CarbohydrateLink:
    @property
    def donor(self) -> int: ...
    @property
    def acceptor(self) -> int | None: ...
    @property
    def donor_atom(self) -> int: ...
    @property
    def acceptor_atom(self) -> int: ...
    @property
    def provenance(self) -> _BondProvenance: ...

@final
class CarbohydrateReport:
    @property
    def monosaccharides(self) -> list[Monosaccharide]: ...
    @property
    def links(self) -> list[CarbohydrateLink]: ...
    @property
    def terminal_links(self) -> list[CarbohydrateLink]: ...
    @property
    def incomplete_residues(self) -> list[int]: ...
    @property
    def dictionary_version(self) -> str | None: ...

def snfg_symbol(component: str) -> SnfgSymbol | None: ...
def carbohydrates(
    structure: Structure,
    components: str | PathLike[str] | None = None,
    *,
    version: str = "unversioned",
    spatial_fallback: bool = True,
) -> CarbohydrateReport: ...
def polymer_atom_roles() -> dict[str, int]:
    """Polymer atom roles by name, as the integer codes a role rule combines with ``|``."""

@final
class Molecule:
    """One MOL block or SDF record: a molecular graph, its header and its data fields."""

    def __init__(
        self,
        elements: Sequence[str],
        coordinates: Float32Array,
        bonds: UInt32Array | None = None,
        *,
        name: str = "",
        formal_charges: Sequence[int] | None = None,
        properties: Mapping[str, str] | None = None,
    ) -> None:
        """Build a molecule from element symbols, coordinates and ``(first, second, order)`` bonds.

        ``order`` is the MDL order: 1, 2, 3, or 4 for aromatic. ``formal_charges`` is one
        integer per atom (0 for none); ``properties`` become SDF data fields.
        """
    @property
    def name(self) -> str: ...
    @property
    def program(self) -> str: ...
    @property
    def comment(self) -> str: ...
    @property
    def version(self) -> Literal["v2000", "v3000"]: ...
    @property
    def atom_count(self) -> int: ...
    @property
    def bond_count(self) -> int: ...
    @property
    def elements(self) -> list[str]: ...
    @property
    def coordinates(self) -> Float32Array: ...
    @property
    def bonds(self) -> UInt32Array:
        """``(first, second, order)`` rows, zero-based; order 4 is aromatic."""
    @property
    def formal_charges(self) -> list[int]: ...
    @property
    def properties(self) -> dict[str, str]: ...
    def to_structure(self) -> Structure:
        """Return the molecule as a structure with its bonds; hydrogens and charges are kept."""
    def to_mol(self) -> str: ...
    def to_sdf(self) -> str: ...

def read_sdf(source: str | bytes | PathLike[str]) -> list[Molecule]:
    """Return every record of an SDF file, with its data fields."""

def read_mol(source: str | bytes | PathLike[str]) -> Molecule:
    """Return one MOL block."""

def write_sdf(molecules: Sequence[Molecule]) -> str:
    """Return SDF text for the molecules, each followed by its data fields."""

def molecule(structure: Structure, *, name: str = "") -> Molecule:
    """Return the graph of a structure: elements, positions and bonds."""

def smarts(structure: Structure, pattern: str) -> list[list[int]]:
    """Return every mapping of a SMARTS pattern onto the structure, as atom indices.

    The structure needs the data the pattern's primitives use (a bond graph; aromaticity,
    formal charges or stereochemistry where the pattern asks), and a missing one is refused.
    """

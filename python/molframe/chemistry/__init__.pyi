from os import PathLike
from typing import Literal, Protocol, final

from .. import Structure

type _RadiusSet = Literal["bondi", "amber_united", "charmm", "alvarez"]

class Float32Array(Protocol):
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

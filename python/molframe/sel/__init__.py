"""Immutable molecular selection expressions."""

from .._native import sel as _native

all = _native.all
aromatic = _native.aromatic
helix = _native.helix
strand = _native.strand
sheet = _native.sheet
alpha_helix = _native.alpha_helix
helix_310 = _native.helix_310
pi_helix = _native.pi_helix
polyproline = _native.polyproline
bridge = _native.bridge
turn = _native.turn
bend = _native.bend
coil = _native.coil
atom = _native.atom
backbone = _native.backbone
chain = _native.chain
glycans = _native.glycans
heavy = _native.heavy
hetero = _native.hetero
hydrogen = _native.hydrogen
ions = _native.ions
ligands = _native.ligands
lipids = _native.lipids
none = _native.none
nonpolar_hydrogen = _native.nonpolar_hydrogen
nucleic = _native.nucleic
nucleic_backbone = _native.nucleic_backbone
nucleic_base = _native.nucleic_base
nucleic_sugar = _native.nucleic_sugar
polar_hydrogen = _native.polar_hydrogen
polymer = _native.polymer
protein = _native.protein
residue = _native.residue
residues_within = _native.residues_within
sidechain = _native.sidechain
water = _native.water
within = _native.within

__all__ = [
    "all",
    "alpha_helix",
    "aromatic",
    "atom",
    "backbone",
    "bend",
    "bridge",
    "chain",
    "coil",
    "glycans",
    "heavy",
    "helix",
    "helix_310",
    "hetero",
    "hydrogen",
    "ions",
    "ligands",
    "lipids",
    "none",
    "nonpolar_hydrogen",
    "nucleic",
    "nucleic_backbone",
    "nucleic_base",
    "nucleic_sugar",
    "pi_helix",
    "polar_hydrogen",
    "polymer",
    "polyproline",
    "protein",
    "residue",
    "residues_within",
    "sheet",
    "sidechain",
    "strand",
    "turn",
    "water",
    "within",
]

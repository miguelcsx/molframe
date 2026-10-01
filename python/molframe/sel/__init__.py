"""Immutable molecular selection expressions."""

from .._native import sel as _native

all = _native.all
aromatic = _native.aromatic
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
nucleic = _native.nucleic
nucleic_backbone = _native.nucleic_backbone
nucleic_base = _native.nucleic_base
nucleic_sugar = _native.nucleic_sugar
polymer = _native.polymer
protein = _native.protein
residue = _native.residue
residues_within = _native.residues_within
sidechain = _native.sidechain
water = _native.water
within = _native.within

__all__ = [
    "all",
    "aromatic",
    "atom",
    "backbone",
    "chain",
    "glycans",
    "heavy",
    "hetero",
    "hydrogen",
    "ions",
    "ligands",
    "lipids",
    "none",
    "nucleic",
    "nucleic_backbone",
    "nucleic_base",
    "nucleic_sugar",
    "polymer",
    "protein",
    "residue",
    "residues_within",
    "sidechain",
    "water",
    "within",
]

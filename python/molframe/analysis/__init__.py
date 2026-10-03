"""Structure analyses."""

from .._native import analysis as _native

ContactTable = _native.ContactTable
GridSpec = _native.GridSpec
ScalarGrid = _native.ScalarGrid
contact_potential = _native.contact_potential
atom_contacts = _native.atom_contacts
contacts = _native.contacts
hydrogen_bonds = _native.hydrogen_bonds
salt_bridges = _native.salt_bridges
dssp = _native.dssp

__all__ = [
    "ContactTable",
    "GridSpec",
    "ScalarGrid",
    "atom_contacts",
    "contact_potential",
    "contacts",
    "dssp",
    "hydrogen_bonds",
    "salt_bridges",
]

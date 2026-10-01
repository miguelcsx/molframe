"""Structure analyses."""
from .._native import analysis as _native

ContactTable = _native.ContactTable
atom_contacts = _native.atom_contacts
contacts = _native.contacts
hydrogen_bonds = _native.hydrogen_bonds
salt_bridges = _native.salt_bridges

__all__ = ["ContactTable", "atom_contacts", "contacts", "hydrogen_bonds", "salt_bridges"]

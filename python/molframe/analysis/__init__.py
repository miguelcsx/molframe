"""Structure analyses."""
from .._native import analysis as _native

ContactTable = _native.ContactTable
atom_contacts = _native.atom_contacts

__all__ = ["ContactTable", "atom_contacts"]

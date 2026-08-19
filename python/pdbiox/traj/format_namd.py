"""NAMD binary coordinate snapshot bindings."""

from .._trajectory import NamdBinary, NamdEndian, NamdError, parse_namd_binary, write_namd_binary

__all__ = ["NamdBinary", "NamdEndian", "NamdError", "parse_namd_binary", "write_namd_binary"]

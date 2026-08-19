"""LAMMPS dump trajectory format."""

from .._trajectory import LammpsError, parse_lammps_dump

__all__ = ["LammpsError", "parse_lammps_dump"]

"""Predicted aligned error and ModelCIF confidence."""

from .._native import modelcif as _native

PaeMatrix = _native.PaeMatrix
pae = _native.pae
read_pae_json = _native.read_pae_json

__all__ = ["PaeMatrix", "pae", "read_pae_json"]

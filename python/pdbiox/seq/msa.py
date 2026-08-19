"""Multiple-sequence-alignment primitives backed by Rust."""

from .._native.seq import MsaError, MsaOptions, multiple_sequence_alignment, progressive_msa

__all__ = ["MsaError", "MsaOptions", "multiple_sequence_alignment", "progressive_msa"]

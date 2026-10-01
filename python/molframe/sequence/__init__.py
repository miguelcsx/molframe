"""Sequence operations."""
from .._native import sequence as _native

Alignment = _native.Alignment
FastaRecord = _native.FastaRecord
align = _native.align
kmer_counts = _native.kmer_counts
parse_fasta = _native.parse_fasta
write_fasta = _native.write_fasta

__all__ = ["Alignment", "FastaRecord", "align", "kmer_counts", "parse_fasta", "write_fasta"]

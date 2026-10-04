"""Sequence operations."""

from .._native import sequence as _native

Alignment = _native.Alignment
FastaRecord = _native.FastaRecord
Scoring = _native.Scoring
SubstitutionMatrix = _native.SubstitutionMatrix
Tree = _native.Tree
align = _native.align
kmer_counts = _native.kmer_counts
msa = _native.msa
neighbor_joining = _native.neighbor_joining
parse_fasta = _native.parse_fasta
upgma = _native.upgma
write_fasta = _native.write_fasta

__all__ = [
    "Alignment",
    "FastaRecord",
    "Scoring",
    "SubstitutionMatrix",
    "Tree",
    "align",
    "kmer_counts",
    "msa",
    "neighbor_joining",
    "parse_fasta",
    "upgma",
    "write_fasta",
]

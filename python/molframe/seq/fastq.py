"""FASTQ records, parsing, and writing."""

from . import FastqError, FastqRecord, parse_fastq, write_fastq

__all__ = ["FastqError", "FastqRecord", "parse_fastq", "write_fastq"]

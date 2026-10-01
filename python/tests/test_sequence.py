"""Pairwise alignment, FASTA and k-mers through the Python contract."""

import pytest

import molframe
from molframe.sequence import FastaRecord


def test_a_global_alignment_pads_the_shorter_sequence_with_gaps():
    result = molframe.sequence.align("GATTACA", "GATACA")
    top, bottom = result.aligned
    assert len(top) == len(bottom) == len(result)
    assert top.replace("-", "") == "GATTACA"
    assert bottom.replace("-", "") == "GATACA"
    # six matches and one single-residue gap opened at -2
    assert result.score == 6 - 2


def test_a_local_alignment_trims_flanking_mismatches():
    result = molframe.sequence.align("TTTTACGTACGTTTTT", "ACGTACGT", mode="local")
    assert result.score == 8
    assert result.aligned == ("ACGTACGT", "ACGTACGT")


def test_unknown_modes_and_positive_gap_costs_are_rejected():
    with pytest.raises(ValueError):
        molframe.sequence.align("A", "A", mode="banded")
    with pytest.raises(ValueError):
        molframe.sequence.align("A", "A", gap_open=1)


def test_fasta_round_trips_and_strips_whitespace():
    text = ">one first record\nACGT\nACGT\n>two\nTTGA\n"
    records = molframe.sequence.parse_fasta(text)
    assert [(r.id, r.description, r.sequence) for r in records] == [
        ("one", "first record", "ACGTACGT"),
        ("two", "", "TTGA"),
    ]
    again = molframe.sequence.parse_fasta(molframe.sequence.write_fasta(records))
    assert [r.sequence for r in again] == [r.sequence for r in records]
    assert len(FastaRecord("x", "ACG")) == 3


def test_kmer_counts_count_repeats_and_reject_zero():
    counts = dict(molframe.sequence.kmer_counts("ACACAC", 2))
    assert counts == {"AC": 3, "CA": 2}
    with pytest.raises(ValueError):
        molframe.sequence.kmer_counts("ACGT", 0)


@pytest.mark.parametrize(
    "name", ["geometry", "surface", "compare", "sequence", "crystal", "analysis", "sel"]
)
def test_every_domain_subpackage_imports_and_matches_the_root_attribute(name):
    import importlib

    package = importlib.import_module(f"molframe.{name}")
    root = getattr(molframe, name)
    for exported in package.__all__:
        assert getattr(package, exported) is getattr(root, exported)

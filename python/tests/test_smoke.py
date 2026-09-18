"""Binding smoke tests (16 §1.1): conversion, identity, and error mapping.

Behaviour lives in Rust and is tested there. These assert only what the
binding layer owns: the module imports fast, arguments convert, borrowed
buffers stay borrowed, and errors carry their findings.
"""

import subprocess
import sys
from pathlib import Path

import numpy
import pytest

import molframe

DATA = Path(__file__).resolve().parents[2] / "crates" / "molframe-py" / "tests" / "data"

# 19 §2: importing the extension stays under 50 ms.
IMPORT_BUDGET_SECONDS = 0.05


def test_import_stays_under_the_budget():
    source = (
        "import time; start = time.perf_counter(); "
        "import molframe; print(time.perf_counter() - start)"
    )
    out = subprocess.run(
        [sys.executable, "-c", source], check=True, capture_output=True, text=True
    ).stdout
    assert float(out) < IMPORT_BUDGET_SECONDS


def test_read_dispatches_by_file_name():
    structure = molframe.read(DATA / "basic.pdb")
    assert structure.atom_count == 2


def test_read_bytes_converts_and_names_the_buffer():
    report = molframe.read_bytes(
        (DATA / "basic.cif").read_bytes(), molframe.ReadOptions.standard(), name="basic.cif"
    )
    assert report.structure.atom_count == 2
    assert report.findings == []


def test_select_accepts_text_and_compiled_query():
    structure = molframe.read(DATA / "basic.pdb")
    by_text = structure.select("name CA")
    by_query = structure.select(molframe.Query("name CA"))
    assert len(by_text) == 1
    assert by_text.indices.tolist() == by_query.indices.tolist()


def test_positions_are_borrowed_not_copied():
    structure = molframe.read(DATA / "basic.pdb")
    first = structure.positions
    second = structure.positions
    assert first.dtype == numpy.float32
    assert numpy.shares_memory(first, second)


def test_errors_carry_their_findings():
    with pytest.raises(molframe.MolframeError) as caught:
        molframe.read_bytes(b"garbage", molframe.ReadOptions.standard(), name="basic.cif")
    error = caught.value
    assert isinstance(error, molframe.SchemaError), "bad CIF syntax must map to a typed error"
    assert error.findings, "the error must carry its findings"
    assert error.code
    assert error.message
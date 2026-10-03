import _thread
import threading
import time

import numpy as np
import pytest

import molframe
from molframe import Cancelled, ExecutionContext, MemoryBudgetError, MolframeError

PDB = """\
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
ATOM      3  C   ALA A   1      13.100   7.600  -6.504  1.00  0.00           C
END
"""


def _cloud(count, seed=1):
    rng = np.random.default_rng(seed)
    return (rng.random((count, 3)) * 40.0).astype(np.float32)


def test_a_context_reports_the_limits_it_was_given():
    context = ExecutionContext(workers=2, memory_budget=1 << 24, scratch_bytes=0)
    assert (context.workers, context.memory_budget, context.scratch_bytes) == (2, 1 << 24, 0)
    assert context.temp_directory is None
    assert not context.is_cancelled
    assert "workers=Some(2)" in repr(context)


def test_inconsistent_limits_are_refused_when_the_context_is_made():
    with pytest.raises(MolframeError):
        ExecutionContext(workers=0)
    with pytest.raises(MolframeError):
        ExecutionContext(memory_budget=0)
    with pytest.raises(MolframeError, match="E6102"):
        ExecutionContext(memory_budget=10, scratch_bytes=1_000)


def test_every_heavy_operation_accepts_a_context_and_agrees_with_the_default():
    coordinates = _cloud(2_500)
    radii = np.full(len(coordinates), 1.7, dtype=np.float32)
    context = ExecutionContext(workers=2)
    default_first, default_second, default_distance = molframe.spatial.neighbor_pairs(
        coordinates, 2.0
    )
    first, second, distance = molframe.spatial.neighbor_pairs(coordinates, 2.0, context=context)
    assert np.array_equal(first, default_first)
    assert np.array_equal(second, default_second)
    assert np.array_equal(distance, default_distance)
    plain = molframe.surface.sasa(coordinates, radii, points=64)
    governed = molframe.surface.sasa(coordinates, radii, points=64, context=context)
    assert np.array_equal(plain, governed)
    structure = molframe.read(PDB.encode(), name="t.pdb")
    assert len(molframe.analysis.atom_contacts(structure, 5.0, context=context)) == 3
    assert molframe.analysis.contacts(structure, 5.0, context=context).status == "complete"


def test_a_memory_budget_the_operation_cannot_meet_is_a_memory_error():
    coordinates = _cloud(4_000)
    tight = ExecutionContext(memory_budget=1_024)
    with pytest.raises(MemoryBudgetError, match="E7001") as raised:
        molframe.spatial.neighbor_pairs(coordinates, 3.0, context=tight)
    assert isinstance(raised.value, MemoryError)
    assert raised.value.code == "MOLFRAME-E7001"


def test_a_cancelled_context_stops_the_operation_with_cancelled():
    coordinates = _cloud(4_000)
    context = ExecutionContext()
    context.cancel()
    assert context.is_cancelled
    with pytest.raises(Cancelled):
        molframe.spatial.neighbor_pairs(coordinates, 3.0, context=context)


def test_cancelling_from_another_thread_reaches_a_running_operation():
    coordinates = _cloud(60_000, seed=3)
    radii = np.full(len(coordinates), 1.7, dtype=np.float32)
    context = ExecutionContext(workers=1)
    timer = threading.Timer(0.2, context.cancel)
    timer.start()
    started = time.perf_counter()
    try:
        with pytest.raises(Cancelled):
            molframe.surface.sasa(coordinates, radii, points=4_000, context=context)
    finally:
        timer.cancel()
    assert time.perf_counter() - started < 20


def test_ctrl_c_interrupts_a_running_operation_as_keyboard_interrupt():
    coordinates = _cloud(60_000, seed=4)
    radii = np.full(len(coordinates), 1.7, dtype=np.float32)
    timer = threading.Timer(0.2, _thread.interrupt_main)
    timer.start()
    started = time.perf_counter()
    try:
        with pytest.raises(KeyboardInterrupt):
            molframe.surface.sasa(
                coordinates, radii, points=4_000, context=ExecutionContext(workers=1)
            )
    finally:
        timer.cancel()
    assert time.perf_counter() - started < 20


def test_the_interpreter_stays_usable_after_an_interrupt():
    context = ExecutionContext()
    coordinates = _cloud(300)
    first, _, _ = molframe.spatial.neighbor_pairs(coordinates, 6.0, context=context)
    again, _, _ = molframe.spatial.neighbor_pairs(coordinates, 6.0, context=context)
    assert np.array_equal(first, again)

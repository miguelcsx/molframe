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


def test_a_memory_refusal_names_the_figures_and_every_knob():
    with pytest.raises(MemoryBudgetError) as raised:
        molframe.spatial.neighbor_pairs(
            _cloud(4_000), 3.0, context=ExecutionContext(memory_budget=1_024)
        )
    message = str(raised.value)
    for needle in (
        "bytes",
        "ExecutionContext(memory_budget",
        "--memory-budget",
        "MOLFRAME_MEMORY_BUDGET",
    ):
        assert needle in message, message


def test_the_default_budget_scales_with_the_machine():
    assert ExecutionContext().memory_budget is None
    assert ExecutionContext().effective_memory_budget >= 100_000_000
    assert ExecutionContext(memory_budget=5_000_000).effective_memory_budget == 5_000_000


def test_distance_matrix_is_guarded_and_offers_f32_condensed_and_blocks():
    coordinates = _cloud(300)
    full = molframe.geometry.distance_matrix(coordinates)
    assert full.shape == (300, 300)
    assert full.dtype == np.float64
    # 300*300*8 = 720 kB does not fit 500 kB; float32 (360 kB) and condensed do.
    tight = ExecutionContext(memory_budget=500_000, scratch_bytes=0)
    with pytest.raises(MemoryBudgetError) as raised:
        molframe.geometry.distance_matrix(coordinates, context=tight)
    message = str(raised.value)
    for needle in (
        "720000 bytes",
        "float32",
        "condensed",
        "distance_matrix_blocks",
        "MOLFRAME_MEMORY_BUDGET",
    ):
        assert needle in message, message
    narrow = molframe.geometry.distance_matrix(coordinates, dtype="float32", context=tight)
    assert narrow.dtype == np.float32
    np.testing.assert_array_equal(narrow, full.astype(np.float32))
    packed = molframe.geometry.distance_matrix(coordinates, condensed=True, context=tight)
    np.testing.assert_array_equal(packed, full[np.triu_indices(300, 1)])
    stacked = [None] * 300
    for start, block in molframe.geometry.distance_matrix_blocks(coordinates, context=tight):
        stacked[start : start + len(block)] = list(block)
    np.testing.assert_array_equal(np.array(stacked), full)
    with pytest.raises(ValueError, match="dtype"):
        molframe.geometry.distance_matrix(coordinates, dtype="float16")

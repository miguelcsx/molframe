"""Batched TM-score and its prefilter: exactness, shapes, errors, GIL release."""

import threading
import time

import numpy as np
import pytest

from molframe import compare

rng = np.random.default_rng(11)


def chain(length: int) -> np.ndarray:
    steps = rng.normal(size=(length, 3))
    steps *= 3.8 / np.linalg.norm(steps, axis=1, keepdims=True)
    return np.cumsum(steps, axis=0).astype(np.float32)


def library() -> list[np.ndarray]:
    out = []
    for _ in range(3):
        base = chain(60)
        out += [base + rng.normal(scale=0.3, size=base.shape).astype(np.float32) for _ in range(3)]
    out += [chain(60) for _ in range(3)]
    return out


def condensed(count: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(count) for j in range(i + 1, count)]


def test_batch_scores_equal_single_pair_scores():
    structures = library()
    scores = compare.tm_scores(structures)
    expected = [
        compare.tm_score(structures[i], structures[j]) for i, j in condensed(len(structures))
    ]
    assert scores.dtype == np.float64
    assert scores.tolist() == expected


def test_explicit_pairs_follow_the_requested_order():
    structures = library()
    pairs = np.array([[3, 1], [0, 8], [5, 2]], dtype=np.uint64)
    scores = compare.tm_scores(structures, pairs=pairs)
    expected = [compare.tm_score(structures[a], structures[b]) for a, b in pairs.tolist()]
    assert scores.tolist() == expected


@pytest.mark.parametrize("threshold", [0.3, 0.6, 0.9])
def test_exact_screen_never_loses_a_pair_at_or_above_the_threshold(threshold):
    structures = library()
    full = compare.tm_scores(structures)
    screen = compare.tm_screen(structures, threshold=threshold)
    kept = full >= threshold
    assert kept.any()
    assert np.array_equal(screen.scores[kept], full[kept])
    assert screen.pruned_coarse == 0
    assert screen.evaluated + screen.pruned_bound == len(full)
    assert int(np.isnan(screen.scores).sum()) == screen.pruned_bound


def test_coarse_stage_reports_its_own_count():
    structures = library()
    screen = compare.tm_screen(structures, threshold=0.5, coarse_ratio=0.5)
    assert screen.evaluated + screen.pruned_bound + screen.pruned_coarse == len(
        condensed(len(structures))
    )
    full = compare.tm_scores(structures)
    kept = full >= 0.5
    assert np.array_equal(screen.scores[kept], full[kept])


def test_invalid_arguments_raise_ergonomic_errors():
    structures = library()
    for bad in (0.0, 1.5, float("nan")):
        with pytest.raises(ValueError, match="finite"):
            compare.tm_screen(structures, threshold=bad)
    with pytest.raises(ValueError, match="finite"):
        compare.tm_screen(structures, threshold=0.5, coarse_ratio=0.0)
    with pytest.raises(ValueError, match="shape"):
        compare.tm_scores(structures, pairs=np.zeros((2, 3), dtype=np.uint64))
    with pytest.raises(ValueError, match="C-contiguous"):
        compare.tm_scores(
            structures,
            pairs=np.asfortranarray(np.array([[0, 1], [0, 2], [1, 2]], dtype=np.uint64)),
        )
    with pytest.raises(ValueError, match="C-contiguous"):
        compare.tm_scores([np.asfortranarray(chain(30)), chain(30)])
    with pytest.raises(ValueError, match="points"):
        compare.tm_scores([chain(30), chain(31)])


def test_the_gil_is_released_while_scoring():
    structures = [chain(300) for _ in range(40)]
    ticks = []
    done = threading.Event()

    def ticker():
        while not done.is_set():
            ticks.append(time.perf_counter())
            time.sleep(0.005)

    thread = threading.Thread(target=ticker)
    thread.start()
    try:
        compare.tm_scores(structures)
    finally:
        done.set()
        thread.join()
    gaps = np.diff(ticks)
    assert len(ticks) > 5
    assert gaps.max() < 0.5

"""Sliced, atom-subset, chunked and streaming trajectory reads.

`data/trajectory/mdanalysis.npz` holds what MDAnalysis decodes from each of the
three fixture files, so every comparison here is exact (`array_equal`).
"""

from pathlib import Path

import numpy as np
import pytest

import molframe
from molframe import trajectory

DATA = Path(__file__).parent / "data" / "trajectory"
FORMATS = ["xtc", "dcd", "trr"]
EXPECTED = np.load(DATA / "mdanalysis.npz")
FRAMES, ATOMS = EXPECTED["xtc"].shape[:2]


def step(fmt):
    """Picoseconds between frames: the DCD writer records its own 1 ps default."""
    return 1.0 if fmt == "dcd" else 2.0


def path(fmt):
    return DATA / f"walk.{fmt}"


@pytest.fixture(params=FORMATS)
def fmt(request):
    return request.param


def test_a_whole_file_read_matches_mdanalysis_bit_for_bit(fmt):
    data = trajectory.read(path(fmt))
    assert data.format == fmt
    assert data.positions.dtype == np.float32
    assert np.array_equal(data.positions, EXPECTED[fmt])
    assert np.array_equal(data.frame_indices, np.arange(FRAMES))
    assert np.allclose(data.times, np.arange(FRAMES) * step(fmt), rtol=1e-6)


@pytest.mark.parametrize(
    "frames",
    [slice(None, None, 3), slice(2, 8), slice(1, None, 4), slice(-4, None), slice(None, -7, 2)],
)
def test_a_frame_slice_keeps_the_original_numbering(fmt, frames):
    data = trajectory.read(path(fmt), frames=frames)
    wanted = np.arange(FRAMES)[frames]
    assert np.array_equal(data.frame_indices, wanted)
    assert np.array_equal(data.positions, EXPECTED[fmt][wanted])
    assert np.allclose(data.times, wanted * step(fmt), rtol=1e-6)


def test_an_atom_subset_follows_the_order_asked(fmt):
    atoms = [5, 0, 17]
    data = trajectory.read(path(fmt), atoms=atoms, frames=slice(None, None, 2))
    assert data.n_atoms == 3
    assert np.array_equal(data.positions, EXPECTED[fmt][::2][:, atoms])


def test_a_slice_past_the_end_is_empty(fmt):
    data = trajectory.read(path(fmt), frames=slice(50, 60))
    assert data.n_frames == 0
    assert data.positions.shape == (0, ATOMS, 3)


def test_a_zero_or_reversed_stride_is_refused(fmt):
    with pytest.raises(ValueError, match="stride"):
        trajectory.read(path(fmt), frames=slice(None, None, 0))
    with pytest.raises(ValueError, match="stride"):
        trajectory.read(path(fmt), frames=slice(None, None, -1))


def test_an_atom_outside_the_frame_is_refused(fmt):
    with pytest.raises(ValueError, match=r"out of range|selected atom"):
        trajectory.read(path(fmt), atoms=[ATOMS])


def test_a_memory_limit_too_small_for_one_frame_is_refused(fmt):
    with pytest.raises(Exception, match=r"limit|storage|memory"):
        trajectory.read(path(fmt), memory_limit=16)


@pytest.mark.parametrize("size", [1, 3, 4, 10, 64])
def test_chunks_concatenate_to_the_whole_read(fmt, size):
    chunks = list(trajectory.read_chunks(path(fmt), chunk_frames=size))
    assert [chunk.n_frames for chunk in chunks[:-1]] == [size] * (len(chunks) - 1)
    assert 0 < chunks[-1].n_frames <= size
    assert np.array_equal(np.concatenate([c.positions for c in chunks]), EXPECTED[fmt])
    assert np.array_equal(np.concatenate([c.frame_indices for c in chunks]), np.arange(FRAMES))


def test_chunks_respect_a_stride_and_an_atom_subset(fmt):
    chunks = trajectory.read_chunks(
        path(fmt), chunk_frames=2, frames=slice(1, None, 3), atoms=[3, 4]
    )
    assert (chunks.n_atoms, chunks.chunk_frames) == (2, 2)
    assert chunks.n_frames in (None, 3)
    got = np.concatenate([chunk.positions for chunk in chunks])
    assert np.array_equal(got, EXPECTED[fmt][1::3][:, [3, 4]])


def test_chunks_are_exhausted_once_and_can_be_rewound(fmt):
    chunks = trajectory.read_chunks(path(fmt), chunk_frames=4)
    assert sum(chunk.n_frames for chunk in chunks) == FRAMES
    assert next(chunks, None) is None
    chunks.seek(6)
    assert np.array_equal(next(chunks).positions, EXPECTED[fmt][6:10])
    chunks.seek(0)
    assert np.array_equal(next(chunks).positions, EXPECTED[fmt][0:4])


def test_an_unrecognised_chunk_size_is_refused():
    with pytest.raises(ValueError, match="chunk_frames"):
        trajectory.read_chunks(path("xtc"), chunk_frames=0)


def test_a_format_without_a_pull_reader_still_reads_sliced(tmp_path):
    text = "".join(f"2\nframe {i}\nC 0 {i} 0\nO 1 {i} 0\n" for i in range(5))
    xyz = tmp_path / "walk.xyz"
    xyz.write_text(text)
    data = trajectory.read(xyz, frames=slice(1, None, 2), atoms=[1])
    assert np.array_equal(data.frame_indices, [1, 3])
    assert np.allclose(data.positions[:, 0, 1], [1, 3])
    with pytest.raises(ValueError, match="pull-based"):
        trajectory.read_chunks(xyz)


def test_the_streamed_rmsd_equals_the_dense_rmsd(fmt):
    dense = trajectory.rmsd(EXPECTED[fmt], reference=0, align=True).value
    series = trajectory.rmsd_stream(path(fmt), 0, align=True)
    assert np.array_equal(series.frames, np.arange(FRAMES))
    assert np.allclose(series.times, np.arange(FRAMES) * step(fmt), rtol=1e-6)
    assert np.allclose(series.values, dense, atol=1e-4)
    by_array = trajectory.rmsd_stream(path(fmt), EXPECTED[fmt][0], align=True)
    assert np.array_equal(by_array.values, series.values)
    assert len(series) == FRAMES


def test_a_streamed_rmsd_honours_a_slice_and_a_subset(fmt):
    series = trajectory.rmsd_stream(
        path(fmt), 2, frames=slice(2, None, 4), atoms=[0, 1, 2, 3, 4, 5], align=False
    )
    positions = EXPECTED[fmt][2::4][:, :6]
    reference = EXPECTED[fmt][2, :6]
    wanted = np.sqrt(((positions - reference) ** 2).sum(-1).mean(-1))
    assert np.array_equal(series.frames, [2, 6])
    assert np.allclose(series.values, wanted, atol=1e-5)


def test_the_streamed_rmsf_matches_numpy(fmt):
    rmsf = trajectory.rmsf_stream(path(fmt))
    positions = EXPECTED[fmt].astype(np.float64)
    wanted = np.sqrt(positions.var(axis=0).sum(-1))
    assert rmsf.shape == (ATOMS,)
    assert np.allclose(rmsf, wanted, atol=1e-9)


def test_streamed_contacts_count_minimum_image_pairs(fmt):
    cell = np.array([30.0, 31.0, 32.0])
    series = trajectory.contact_counts_stream(path(fmt), 6.0)
    assert series.values.dtype == np.uint64
    for frame, count in zip(series.frames, series.values, strict=True):
        positions = EXPECTED[fmt][frame].astype(np.float64)
        delta = positions[:, None] - positions[None]
        delta -= cell * np.round(delta / cell)
        close = (delta**2).sum(-1) <= 36.0
        assert count == (close.sum() - ATOMS) // 2


def test_a_callback_stream_sees_each_frame_once_and_stops_on_an_exception(fmt):
    seen = []

    def collect(frame, time, positions):
        assert not positions.flags.writeable
        seen.append((frame, time, positions))

    count = trajectory.run_analysis_stream(path(fmt), collect, frames=slice(None, None, 3))
    assert count == len(seen) == 4
    assert [frame for frame, _, _ in seen] == [0, 3, 6, 9]
    assert all(np.array_equal(positions, EXPECTED[fmt][frame]) for frame, _, positions in seen)

    def fail(*_: object) -> None:
        message = "stop at two"
        raise RuntimeError(message)

    with pytest.raises(RuntimeError, match="stop at two"):
        trajectory.run_analysis_stream(path(fmt), fail)


def test_a_stream_under_a_cancelled_context_does_not_run():
    context = molframe.ExecutionContext()
    context.cancel()
    with pytest.raises(Exception, match="cancel"):
        trajectory.rmsf_stream(path("xtc"), context=context)


def test_a_stream_budget_too_small_for_the_frame_workspace_is_refused():
    context = molframe.ExecutionContext(memory_budget=1024)
    with pytest.raises(molframe.MemoryBudgetError):
        trajectory.rmsf_stream(path("dcd"), context=context)


def test_float64_and_strided_coordinates_are_accepted_like_the_contiguous_float32():
    positions = np.ascontiguousarray(EXPECTED["xtc"])
    wide = positions.astype(np.float64)
    interleaved = np.repeat(positions, 2, axis=1)[:, ::2]
    assert not interleaved.flags.c_contiguous
    transposed = np.asfortranarray(positions)
    base = trajectory.pairwise_rmsd(positions).value
    for variant in (wide, interleaved, transposed):
        assert np.array_equal(trajectory.pairwise_rmsd(variant).value, base)
        assert np.array_equal(trajectory.rmsd(variant).value, trajectory.rmsd(positions).value)


def test_coordinates_that_are_not_frames_by_atoms_by_three_are_refused():
    with pytest.raises(ValueError, match="shape"):
        trajectory.rmsd(np.zeros((4, 5), dtype=np.float32))
    with pytest.raises(ValueError, match="shape"):
        trajectory.pairwise_rmsd(np.zeros((4, 5, 2), dtype=np.float32))


def blocks(fmt, workers, size=3):
    context = molframe.ExecutionContext(workers=workers)
    return list(trajectory.read_chunks(path(fmt), chunk_frames=size, context=context))


@pytest.mark.parametrize("workers", [1, 2, 4, 8])
def test_parallel_blocks_equal_the_serial_blocks_in_order(fmt, workers):
    serial = list(trajectory.read_chunks(path(fmt), chunk_frames=3))
    parallel = blocks(fmt, workers)
    assert [c.n_frames for c in parallel] == [c.n_frames for c in serial]
    for got, want in zip(parallel, serial, strict=True):
        assert np.array_equal(got.positions, want.positions)
        assert np.array_equal(got.frame_indices, want.frame_indices)
    assert np.array_equal(np.concatenate([c.frame_indices for c in parallel]), np.arange(FRAMES))


def test_parallel_blocks_honour_a_stride_and_an_atom_subset():
    context = molframe.ExecutionContext(workers=4)
    chunks = trajectory.read_chunks(
        path("xtc"), chunk_frames=2, frames=slice(1, None, 3), atoms=[3, 4], context=context
    )
    got = np.concatenate([chunk.positions for chunk in chunks])
    assert np.array_equal(got, EXPECTED["xtc"][1::3][:, [3, 4]])


def test_a_cancelled_context_stops_the_parallel_stream():
    context = molframe.ExecutionContext(workers=4)
    chunks = trajectory.read_chunks(path("xtc"), chunk_frames=1, context=context)
    next(chunks)
    context.cancel()
    with pytest.raises(molframe.Cancelled):
        list(chunks)

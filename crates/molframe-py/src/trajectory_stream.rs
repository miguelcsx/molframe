//! Frame-sliced, atom-subset and chunked reading of trajectories.
//!
//! XTC, DCD and TRR are pulled frame by frame into one exactly sized output
//! buffer that `NumPy` then owns, so a read peaks at the output plus one frame
//! rather than at two copies of the file. Other containers are materialised by
//! the Rust reader and then sliced.

use crate::trajectory::PyTrajectory;
use molframe::trajectory::{
    DenseChunk, FrameSlice, ParallelChunks, RandomAccess, SelectedReader, Timestep,
    TrajectoryError, TrajectoryFormat, TrajectoryReadOptions, TrajectoryReader,
    TrajectoryReaderOptions, count_frames, read_dense_chunk, read_trajectory,
    read_trajectory_materialized,
};
use pyo3::prelude::*;
use pyo3::types::PySlice;
use std::path::PathBuf;
use std::sync::Mutex;

/// A pull reader restricted to a frame slice and an atom subset.
pub(crate) type Reader = SelectedReader<Box<dyn TrajectoryReader>>;

/// Frame-workspace bytes per source atom: source and output positions, decoder
/// scratch and the read-ahead index, with headroom.
const WORKSPACE_PER_ATOM: usize = 64;

/// Fixed part of the frame workspace: input buffer and bookkeeping.
const WORKSPACE_FIXED: usize = 16 << 20;

/// Where a trajectory lives and how its reader is bounded.
#[derive(Debug)]
pub(crate) struct Source {
    path: PathBuf,
    format: Option<TrajectoryFormat>,
    memory_limit: usize,
}

impl Source {
    pub(crate) fn new(
        path: PathBuf,
        format: Option<&str>,
        memory_limit: Option<usize>,
    ) -> PyResult<Self> {
        let format = format
            .map(|name| name.parse().map_err(crate::error::kernel))
            .transpose()?
            .or_else(|| TrajectoryFormat::infer(&path));
        Ok(Self {
            path,
            format,
            memory_limit: match memory_limit {
                Some(limit) => limit,
                None => TrajectoryReaderOptions::DEFAULT_MEMORY_LIMIT_BYTES,
            },
        })
    }

    fn pulls(&self) -> bool {
        matches!(
            self.format,
            Some(TrajectoryFormat::Xtc | TrajectoryFormat::Dcd | TrajectoryFormat::Trr)
        )
    }

    fn open_raw(&self) -> PyResult<Box<dyn TrajectoryReader>> {
        let mut raw = read_trajectory(
            &self.path,
            &TrajectoryReaderOptions {
                format: self.format,
                memory_limit_bytes: self.memory_limit,
            },
        )
        .map_err(crate::error::failure)?;
        // XTC learns its frame count only once indexed; seeking to the start builds
        // the index so the count (and an exact output reservation) is known.
        if raw.n_frames().is_none() && raw.random_access() != RandomAccess::None {
            raw.seek(0).map_err(crate::error::kernel)?;
        }
        Ok(raw)
    }

    /// Source frames, counted by a decode pass when the container does not know.
    fn total(&self) -> PyResult<usize> {
        let mut raw = self.open_raw()?;
        match raw.n_frames() {
            Some(total) => Ok(total),
            None => count_frames(&mut raw, &mut Timestep::default()).map_err(crate::error::kernel),
        }
    }

    /// Atoms per frame of the trajectory, as its container records them.
    pub(crate) fn atom_count(&self) -> PyResult<usize> {
        Ok(self.open_raw()?.n_atoms())
    }

    /// The selection reader and the frame workspace the streaming analyses reserve.
    pub(crate) fn open(
        &self,
        frames: Option<&Bound<'_, PySlice>>,
        atoms: Option<Vec<usize>>,
    ) -> PyResult<(Reader, usize)> {
        let slice = resolve_slice(frames, || self.total())?;
        self.open_slice(slice, atoms)
    }

    /// As [`Self::open`] for a slice already resolved, so it needs no Python state.
    pub(crate) fn open_slice(
        &self,
        slice: FrameSlice,
        atoms: Option<Vec<usize>>,
    ) -> PyResult<(Reader, usize)> {
        let raw = self.open_raw()?;
        let workspace = raw
            .n_atoms()
            .saturating_mul(WORKSPACE_PER_ATOM)
            .saturating_add(WORKSPACE_FIXED);
        let reader = SelectedReader::new(raw, slice, atoms).map_err(crate::error::kernel)?;
        Ok((reader, workspace))
    }
}

/// A Python `slice` as a forward frame slice; a negative bound counts from the end.
pub(crate) fn resolve_slice(
    slice: Option<&Bound<'_, PySlice>>,
    total: impl FnOnce() -> PyResult<usize>,
) -> PyResult<FrameSlice> {
    let Some(slice) = slice else {
        return Ok(FrameSlice::ALL);
    };
    let part = |name: &str| -> PyResult<Option<isize>> { slice.getattr(name)?.extract() };
    let (start, stop, stride) = (part("start")?, part("stop")?, part("step")?);
    let stride = match stride {
        Some(stride) => stride,
        None => 1,
    };
    if stride < 1 {
        return Err(crate::error::value(
            "frame stride must be a positive integer; reverse order is not supported",
        ));
    }
    let counted = if start.is_some_and(|bound| bound < 0) || stop.is_some_and(|bound| bound < 0) {
        Some(total()?)
    } else {
        None
    };
    let place = |bound: isize| {
        if let Ok(index) = usize::try_from(bound) {
            return index;
        }
        let known = match counted {
            Some(known) => known,
            None => 0,
        };
        known.saturating_sub(bound.unsigned_abs())
    };
    Ok(FrameSlice {
        start: start.map_or(0, place),
        stop: stop.map_or(usize::MAX, place),
        step: match usize::try_from(stride) {
            Ok(every) => every,
            Err(_) => 1,
        },
    })
}

/// Reads a trajectory, or a slice and atom subset of it, into arrays.
///
/// `frames` is a `slice` of source frames (`slice(0, None, 10)` keeps every
/// tenth); `atoms` lists the atom indices to keep, in the order given. XTC seeks
/// to each requested frame; DCD and TRR decode and discard the frames between.
/// `memory_limit` bounds the reader's own workspace in bytes. The result's
/// `frame_indices` name the source frame of every row.
///
/// A `.pdb`, `.cif`, `.mmcif` or `.bcif` path (or the matching `format`) reads
/// the structure and returns its models as frames; a ragged ensemble is refused
/// with an error naming `ragged_models`.
#[pyfunction]
#[pyo3(signature = (path, *, format=None, frames=None, atoms=None, memory_limit=None))]
pub(crate) fn read(
    py: Python<'_>,
    path: PathBuf,
    format: Option<&str>,
    frames: Option<&Bound<'_, PySlice>>,
    atoms: Option<Vec<usize>>,
    memory_limit: Option<usize>,
) -> PyResult<PyTrajectory> {
    if let Some(name) = crate::trajectory_structure::structure_format(&path, format) {
        return crate::trajectory_structure::read_models(py, &path, name, frames, atoms.as_deref());
    }
    let source = Source::new(path, format, memory_limit)?;
    if let (true, Some(format)) = (source.pulls(), source.format) {
        let (mut reader, _) = source.open(frames, atoms)?;
        // A container that cannot state its length (TRR) is counted by one decode
        // pass so the output is reserved exactly instead of growing by doubling.
        let wanted = if reader.n_frames().is_none() {
            let slice = resolve_slice(frames, || source.total())?;
            slice.count(source.total()?)
        } else {
            usize::MAX
        };
        let chunk = py
            .detach(move || read_dense_chunk(&mut reader, &mut Timestep::default(), wanted))
            .map_err(crate::error::kernel)?;
        return PyTrajectory::from_chunk(py, format.name().replace('-', "_"), chunk);
    }
    let options = TrajectoryReadOptions {
        format: source.format,
        ..TrajectoryReadOptions::default()
    };
    let data = py
        .detach(|| read_trajectory_materialized(&source.path, &options))
        .map_err(crate::error::failure)?;
    let slice = resolve_slice(frames, || Ok(data.frames.len()))?;
    let source_atoms = data.frames.first().map_or(0, |frame| frame.positions.len());
    if let Some(index) = atoms.iter().flatten().find(|index| **index >= source_atoms) {
        return Err(crate::error::index(format!(
            "atom index {index} is out of range for {source_atoms} atoms"
        )));
    }
    let mut chunk = DenseChunk {
        atoms: atoms.as_ref().map_or(source_atoms, Vec::len),
        ..DenseChunk::default()
    };
    let stop = slice.stop.min(data.frames.len());
    for index in (slice.start..stop).step_by(slice.step) {
        let Some(frame) = data.frames.get(index) else {
            break;
        };
        let mut selected = Timestep {
            frame: index,
            time: frame.time,
            cell: frame.cell,
            ..Timestep::default()
        };
        match &atoms {
            Some(atoms) => selected
                .positions
                .extend(atoms.iter().filter_map(|atom| frame.positions.get(*atom))),
            None => selected.positions.extend_from_slice(&frame.positions),
        }
        chunk.push(&selected).map_err(crate::error::kernel)?;
    }
    PyTrajectory::from_chunk(py, data.format.name().replace('-', "_"), chunk)
}

/// How the next block is produced: pulled here, or by background workers.
enum Pull {
    Serial(Box<(Reader, Timestep)>),
    Parallel(ParallelChunks),
}

/// Successive blocks of frames read from one trajectory, each a `Trajectory`.
#[pyclass(
    name = "TrajectoryChunks",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
struct PyTrajectoryChunks {
    state: Mutex<Pull>,
    format: String,
    chunk_frames: usize,
    atoms: usize,
    frames: Option<usize>,
}

impl std::fmt::Debug for PyTrajectoryChunks {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("PyTrajectoryChunks")
            .field("format", &self.format)
            .field("chunk_frames", &self.chunk_frames)
            .field("atoms", &self.atoms)
            .field("frames", &self.frames)
            .finish_non_exhaustive()
    }
}

#[pymethods]
impl PyTrajectoryChunks {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self, py: Python<'_>) -> PyResult<Option<PyTrajectory>> {
        let chunk = py.detach(|| {
            let Ok(mut state) = self.state.lock() else {
                return Err(None);
            };
            match &mut *state {
                Pull::Serial(pair) => {
                    let (reader, frame) = &mut **pair;
                    read_dense_chunk(reader, frame, self.chunk_frames).map_err(Some)
                }
                Pull::Parallel(chunks) => match chunks.next() {
                    Some(chunk) => chunk.map_err(Some),
                    None => Ok(DenseChunk::default()),
                },
            }
        });
        let chunk = match chunk {
            Ok(chunk) => chunk,
            Err(Some(error)) => return Err(crate::error::kernel(error)),
            Err(None) => return Err(crate::error::internal("the chunk reader was poisoned")),
        };
        if chunk.frames() == 0 {
            return Ok(None);
        }
        PyTrajectory::from_chunk(py, self.format.clone(), chunk).map(Some)
    }

    /// Moves to the `frame`-th selected frame; DCD and TRR can only move forward.
    fn seek(&self, frame: usize) -> PyResult<()> {
        let mut state = self
            .state
            .lock()
            .map_err(|_| crate::error::internal("the chunk reader was poisoned"))?;
        match &mut *state {
            Pull::Serial(pair) => pair.0.seek(frame).map_err(crate::error::kernel),
            Pull::Parallel(_) => Err(crate::error::value(
                "seek is not available on a chunk iterator that decodes in parallel",
            )),
        }
    }

    /// Selected frames in the whole iteration, or `None` when the container cannot say.
    #[getter]
    const fn n_frames(&self) -> Option<usize> {
        self.frames
    }

    /// Atoms in every chunk.
    #[getter]
    const fn n_atoms(&self) -> usize {
        self.atoms
    }

    /// Frames per chunk; the last chunk may hold fewer.
    #[getter]
    const fn chunk_frames(&self) -> usize {
        self.chunk_frames
    }
}

/// Iterates a trajectory in blocks of `chunk_frames` frames.
///
/// Memory stays at one block plus one frame however long the trajectory is.
/// `frames`, `atoms` and `memory_limit` mean what they do for `read`. With a
/// `context` whose worker budget exceeds one, background workers decode
/// blocks ahead of the consumer: blocks still arrive in frame order and are
/// identical at every worker count, the read-ahead window is reserved from the
/// context's memory budget, and the context's `cancel()` stops the workers.
#[pyfunction]
#[pyo3(signature = (path, *, chunk_frames=256, format=None, frames=None, atoms=None, memory_limit=None, context=None))]
fn read_chunks(
    py: Python<'_>,
    path: PathBuf,
    chunk_frames: usize,
    format: Option<&str>,
    frames: Option<&Bound<'_, PySlice>>,
    atoms: Option<Vec<usize>>,
    memory_limit: Option<usize>,
    context: Option<&crate::execution::PyExecutionContext>,
) -> PyResult<PyTrajectoryChunks> {
    if chunk_frames == 0 {
        return Err(crate::error::value("chunk_frames must be positive"));
    }
    let source = Source::new(path, format, memory_limit)?;
    let slice = resolve_slice(frames, || source.total())?;
    let (reader, _) = py.detach(|| source.open_slice(slice, atoms.clone()))?;
    let (atom_count, frame_count) = (reader.n_atoms(), reader.n_frames());
    let format = source
        .format
        .map_or_else(String::new, |format| format.name().replace('-', "_"));
    let parallel = match (context, frame_count) {
        (Some(context), Some(total)) if reader.random_access() != RandomAccess::None => {
            let built = context.build(context.token())?;
            (built.worker_budget() > 1 && total > 0).then_some((built, total))
        }
        _ => None,
    };
    let state = match parallel {
        Some((built, total)) => {
            drop(reader);
            let opened = py.detach(move || {
                open_parallel(
                    &source,
                    slice,
                    atoms.as_deref(),
                    total,
                    chunk_frames,
                    &built,
                )
            })?;
            Pull::Parallel(opened)
        }
        None => Pull::Serial(Box::new((reader, Timestep::default()))),
    };
    Ok(PyTrajectoryChunks {
        format,
        chunk_frames,
        atoms: atom_count,
        frames: frame_count,
        state: Mutex::new(state),
    })
}

/// Starts the workers; a reader that fails to open reports its Python error.
fn open_parallel(
    source: &Source,
    slice: FrameSlice,
    atoms: Option<&[usize]>,
    total: usize,
    chunk_frames: usize,
    context: &molframe::ExecutionContext,
) -> PyResult<ParallelChunks> {
    let failure = std::cell::RefCell::new(None);
    let opened = ParallelChunks::open(
        || match source.open_slice(slice, atoms.map(<[usize]>::to_vec)) {
            Ok((reader, _)) => Ok(reader),
            Err(error) => {
                *failure.borrow_mut() = Some(error);
                Err(TrajectoryError::RandomAccessUnavailable)
            }
        },
        0..total,
        chunk_frames,
        context,
    );
    match (opened, failure.take()) {
        (Ok(chunks), _) => Ok(chunks),
        (Err(_), Some(error)) => Err(error),
        (Err(error), None) => Err(crate::error::kernel(error)),
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyTrajectoryChunks>()?;
    module.add_function(wrap_pyfunction!(read_chunks, module)?)?;
    crate::trajectory_analyses::register(module)
}

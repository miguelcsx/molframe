//! Reusable-buffer trajectory reader contract and concrete sources.

use crate::{Frame, Trajectory, TrajectoryBuildError};
use molframe_core::structure::UnitCell;
use std::collections::BTreeMap;

/// Native units reported by a coordinate source before boundary conversion.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Units {
    /// Position unit name.
    pub length: &'static str,
    /// Time unit name.
    pub time: &'static str,
    /// Force unit name.
    pub force: &'static str,
}

impl Units {
    /// molframe canonical units: ångström, picosecond and kJ mol⁻¹ Å⁻¹.
    pub const CANONICAL: Self = Self {
        length: "angstrom",
        time: "picosecond",
        force: "kilojoule_per_mole_angstrom",
    };
}

/// Seeking support promised by a reader.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RandomAccess {
    /// Every frame offset is directly known.
    Full,
    /// Seeking becomes available after building an offset index.
    ViaIndex,
    /// This source can only move forward.
    None,
}

/// Auxiliary per-frame value.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum FrameValue {
    /// Floating-point scalar.
    Float(f64),
    /// Signed integer scalar.
    Integer(i64),
    /// Textual metadata.
    Text(Box<str>),
    /// Floating-point vector.
    Floats(Vec<f64>),
}

/// Reusable coordinate and auxiliary buffers for one frame.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Timestep {
    /// Zero-based frame ordinal.
    pub frame: usize,
    /// Simulation time in picoseconds.
    pub time: Option<f64>,
    /// Time step in picoseconds.
    pub dt: Option<f64>,
    /// Positions in ångström.
    pub positions: Vec<[f32; 3]>,
    /// Velocities in ångström per picosecond.
    pub velocities: Option<Vec<[f32; 3]>>,
    /// Forces in kJ mol⁻¹ Å⁻¹.
    pub forces: Option<Vec<[f32; 3]>>,
    /// Periodic unit cell.
    pub cell: Option<UnitCell>,
    /// Format-specific auxiliary streams in lexical key order.
    pub data: BTreeMap<Box<str>, FrameValue>,
}

impl Timestep {
    fn replace_positions(&mut self, frame: usize, positions: &[[f32; 3]]) {
        self.frame = frame;
        self.positions.clear();
        self.positions.extend_from_slice(positions);
        self.time = None;
        self.dt = None;
        self.velocities = None;
        self.forces = None;
        self.cell = None;
        self.data.clear();
    }
}

/// Reader or topology-pairing failure.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum TrajectoryError {
    /// An owned fixed-width trajectory could not be constructed.
    #[error(transparent)]
    Build(#[from] TrajectoryBuildError),
    /// Seeking was requested from a forward-only source.
    #[error("trajectory source does not support random access")]
    RandomAccessUnavailable,
    /// Sources in a chain disagree on atom count.
    #[error("trajectory atom-count mismatch: expected {expected}, found {found}")]
    AtomCountMismatch {
        /// Atom count established by the first source.
        expected: usize,
        /// Atom count reported by a later source.
        found: usize,
    },
    /// A transform selected an atom outside the frame.
    #[error("trajectory transform selected atom {index}, but the frame has {atoms} atoms")]
    SelectionOutOfRange {
        /// Invalid zero-based atom index.
        index: usize,
        /// Current frame atom count.
        atoms: usize,
    },
    /// A rigid fit had fewer than three non-collinear points.
    #[error("trajectory frame does not define the requested rigid fit")]
    DegenerateFit,
    /// Parallel execution was requested for an order-dependent analysis.
    #[error("frame analysis does not opt into parallel execution")]
    AnalysisNotParallel,
    /// A worker count of zero was requested.
    #[error("frame analysis requires at least one worker")]
    InvalidWorkerCount,
    /// A scoped worker panicked before returning its block results.
    #[error("frame-analysis worker panicked")]
    WorkerPanicked,
    /// The monotonic coordinate-generation identifier reached its representable limit.
    #[error("trajectory coordinate generation is exhausted")]
    CoordinateGenerationExhausted,
    /// A transform produced a coordinate outside finite single precision.
    #[error("trajectory transform produced an unrepresentable coordinate")]
    UnrepresentableCoordinate,
    /// An in-memory trajectory failed to provide an index below its declared length.
    #[error("trajectory omitted declared frame {index}")]
    MissingFrame {
        /// Zero-based frame index requested by the executor.
        index: usize,
    },
    /// The shared spatial backend rejected a cutoff, atom or periodic cell.
    #[error("trajectory spatial operation failed: {0}")]
    Spatial(#[from] molframe_spatial::SpatialError),
    /// A periodic transform was requested for a frame without a unit cell.
    #[error("trajectory frame has no periodic unit cell")]
    MissingCell,
    /// One atom occurred in more than one wrap group.
    #[error("atom {index} occurs in more than one periodic wrap group")]
    OverlappingGroups {
        /// Duplicated atom index.
        index: usize,
    },
    /// A format reader encountered an operating-system I/O failure.
    #[error("{format} trajectory I/O failed: {kind:?}")]
    SourceIo {
        /// Stable source format name.
        format: &'static str,
        /// Portable I/O failure class.
        kind: std::io::ErrorKind,
    },
    /// A format reader rejected a malformed or contradictory frame.
    #[error("invalid {format} trajectory record")]
    InvalidSource {
        /// Stable source format name.
        format: &'static str,
    },
    /// Reader-owned and caller-output storage would exceed the explicit ceiling.
    #[error("trajectory frame storage requires {required} bytes, over the {limit} byte limit")]
    MemoryLimit {
        /// Required bytes, or `usize::MAX` when the dimension overflows.
        required: usize,
        /// Caller-provided ceiling.
        limit: usize,
    },
    /// Cooperative cancellation stopped the pull before consuming another frame.
    #[error("trajectory batch source was cancelled")]
    Cancelled,
    /// A global frame or chunk identity reached the end of its 64-bit domain.
    #[error("trajectory batch identity overflow")]
    IdentityOverflow,
}

/// Coordinate source that fills a caller-owned timestep buffer.
pub trait TrajectoryReader: Send {
    /// Stable format or source name.
    fn format(&self) -> &'static str;

    /// Number of atoms in every frame.
    fn n_atoms(&self) -> usize;
    /// Frame count, absent for open-ended streams.
    fn n_frames(&self) -> Option<usize>;
    /// Source units after reader-boundary conversion.
    fn units(&self) -> Units;
    /// Declared seeking capability.
    fn random_access(&self) -> RandomAccess;
    /// Reads the next frame into a reusable buffer; `false` means EOF.
    ///
    /// # Errors
    ///
    /// Returns a source or atom-count error without publishing a partial frame.
    fn read_next(&mut self, timestep: &mut Timestep) -> Result<bool, TrajectoryError>;
    /// Reads after validating the complete frame workspace against a byte ceiling.
    ///
    /// Readers without a bounded decoder refuse before advancing. The ceiling
    /// includes retained output capacity and decoder workspace, not the source
    /// object supplied by the caller.
    ///
    /// # Errors
    ///
    /// Returns a resource or capability error before an unsupported read.
    fn read_next_bounded(
        &mut self,
        _timestep: &mut Timestep,
        _bytes: usize,
    ) -> Result<bool, TrajectoryError> {
        Err(TrajectoryError::InvalidSource {
            format: "reader has no bounded decoder",
        })
    }

    /// Positions the next read at one frame.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryError::RandomAccessUnavailable`] for streams.
    fn seek(&mut self, frame: usize) -> Result<(), TrajectoryError>;
    /// Positions the next read at the first frame.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryError::RandomAccessUnavailable`] for streams.
    fn rewind(&mut self) -> Result<(), TrajectoryError> {
        self.seek(0)
    }
}

impl<R: TrajectoryReader + ?Sized> TrajectoryReader for Box<R> {
    fn format(&self) -> &'static str {
        (**self).format()
    }

    fn n_atoms(&self) -> usize {
        (**self).n_atoms()
    }

    fn n_frames(&self) -> Option<usize> {
        (**self).n_frames()
    }

    fn units(&self) -> Units {
        (**self).units()
    }

    fn random_access(&self) -> RandomAccess {
        (**self).random_access()
    }

    fn read_next(&mut self, timestep: &mut Timestep) -> Result<bool, TrajectoryError> {
        (**self).read_next(timestep)
    }

    fn read_next_bounded(
        &mut self,
        timestep: &mut Timestep,
        bytes: usize,
    ) -> Result<bool, TrajectoryError> {
        (**self).read_next_bounded(timestep, bytes)
    }

    fn seek(&mut self, frame: usize) -> Result<(), TrajectoryError> {
        (**self).seek(frame)
    }
}

/// Random-access reader over an explicitly materialised trajectory.
#[derive(Debug)]
pub struct MemoryReader<'a> {
    trajectory: &'a Trajectory,
    next: usize,
    atoms: usize,
}

impl<'a> MemoryReader<'a> {
    /// Wraps a memory trajectory without copying its frames.
    #[must_use]
    pub fn new(trajectory: &'a Trajectory) -> Self {
        Self {
            trajectory,
            next: 0,
            atoms: trajectory.frame(0).map_or(0, |frame| frame.positions.len()),
        }
    }
}

impl TrajectoryReader for MemoryReader<'_> {
    fn format(&self) -> &'static str {
        "memory"
    }

    fn n_atoms(&self) -> usize {
        self.atoms
    }

    fn n_frames(&self) -> Option<usize> {
        Some(self.trajectory.len())
    }

    fn units(&self) -> Units {
        Units::CANONICAL
    }

    fn random_access(&self) -> RandomAccess {
        RandomAccess::Full
    }

    fn read_next(&mut self, timestep: &mut Timestep) -> Result<bool, TrajectoryError> {
        let Some(frame) = self.trajectory.frame(self.next) else {
            return Ok(false);
        };
        timestep.replace_positions(self.next, frame.positions);
        self.next += 1;
        Ok(true)
    }

    fn read_next_bounded(
        &mut self,
        timestep: &mut Timestep,
        bytes: usize,
    ) -> Result<bool, TrajectoryError> {
        super::stream_consume::prepare_positions(timestep, self.n_atoms(), bytes)?;
        self.read_next(timestep)
    }

    fn seek(&mut self, frame: usize) -> Result<(), TrajectoryError> {
        self.next = frame.min(self.trajectory.len());
        Ok(())
    }
}

/// Forward-only source that owns only its iterator and the current frame.
#[derive(Debug)]
pub struct StreamingReader<I> {
    frames: I,
    atoms: usize,
    next: usize,
}

impl<I> StreamingReader<I> {
    /// Creates a stream with an explicit fixed atom count.
    #[must_use]
    pub const fn new(frames: I, atoms: usize) -> Self {
        Self {
            frames,
            atoms,
            next: 0,
        }
    }
}

impl<I: Iterator<Item = Frame> + Send> TrajectoryReader for StreamingReader<I> {
    fn format(&self) -> &'static str {
        "stream"
    }

    fn n_atoms(&self) -> usize {
        self.atoms
    }

    fn n_frames(&self) -> Option<usize> {
        None
    }

    fn units(&self) -> Units {
        Units::CANONICAL
    }

    fn random_access(&self) -> RandomAccess {
        RandomAccess::None
    }

    fn read_next(&mut self, timestep: &mut Timestep) -> Result<bool, TrajectoryError> {
        let Some(frame) = self.frames.next() else {
            return Ok(false);
        };
        if frame.positions.len() != self.atoms {
            return Err(TrajectoryError::AtomCountMismatch {
                expected: self.atoms,
                found: frame.positions.len(),
            });
        }
        timestep.replace_positions(self.next, &frame.positions);
        self.next += 1;
        Ok(true)
    }

    fn read_next_bounded(
        &mut self,
        timestep: &mut Timestep,
        bytes: usize,
    ) -> Result<bool, TrajectoryError> {
        super::stream_consume::prepare_positions(timestep, self.n_atoms(), bytes)?;
        self.read_next(timestep)
    }

    fn seek(&mut self, _frame: usize) -> Result<(), TrajectoryError> {
        Err(TrajectoryError::RandomAccessUnavailable)
    }
}

/// Several homogeneous readers presented as one forward sequence.
#[derive(Debug)]
pub struct ChainedReader<R> {
    readers: Vec<R>,
    current: usize,
    frame: usize,
    atoms: usize,
}

impl<R: TrajectoryReader> ChainedReader<R> {
    /// Validates atom counts before consuming any source.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryError::AtomCountMismatch`] if sources cannot share a
    /// topology.
    pub fn new(readers: Vec<R>) -> Result<Self, TrajectoryError> {
        let atoms = readers.first().map_or(0, TrajectoryReader::n_atoms);
        for reader in &readers {
            if reader.n_atoms() != atoms {
                return Err(TrajectoryError::AtomCountMismatch {
                    expected: atoms,
                    found: reader.n_atoms(),
                });
            }
        }
        Ok(Self {
            readers,
            current: 0,
            frame: 0,
            atoms,
        })
    }
}

impl<R: TrajectoryReader> TrajectoryReader for ChainedReader<R> {
    fn format(&self) -> &'static str {
        "chain"
    }

    fn n_atoms(&self) -> usize {
        self.atoms
    }

    fn n_frames(&self) -> Option<usize> {
        self.readers.iter().try_fold(0usize, |total, reader| {
            reader.n_frames().and_then(|count| total.checked_add(count))
        })
    }

    fn units(&self) -> Units {
        Units::CANONICAL
    }

    fn random_access(&self) -> RandomAccess {
        RandomAccess::None
    }

    fn read_next(&mut self, timestep: &mut Timestep) -> Result<bool, TrajectoryError> {
        while let Some(reader) = self.readers.get_mut(self.current) {
            if reader.read_next(timestep)? {
                timestep.frame = self.frame;
                self.frame += 1;
                return Ok(true);
            }
            self.current += 1;
        }
        Ok(false)
    }

    fn seek(&mut self, _frame: usize) -> Result<(), TrajectoryError> {
        Err(TrajectoryError::RandomAccessUnavailable)
    }
}

#[cfg(test)]
#[path = "reader_tests.rs"]
mod tests;

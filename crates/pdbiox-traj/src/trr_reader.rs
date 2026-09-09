//! Forward TRR decoding with reusable output arrays and bounded buffering.

use super::{Header, TrrPrecision, Xdr, read_header};
use crate::{RandomAccess, Timestep, TrajectoryError, TrajectoryReader, Units};
use std::fs::File;
use std::io::{BufReader, Read, Seek};
use std::mem::size_of;
use std::path::Path;

const BUFFER: usize = 8192;

/// A forward-only TRR reader preserving positions, velocities, forces and cell.
#[derive(Debug)]
pub struct TrrReader {
    input: BufReader<File>,
    pending: Option<Header>,
    atoms: usize,
    next: usize,
    previous_time: Option<f64>,
    last_step: Option<i32>,
    last_precision: Option<TrrPrecision>,
    limit: usize,
    file_len: u64,
}

impl TrrReader {
    /// Opens a TRR source without loading its frames.
    ///
    /// # Errors
    ///
    /// Returns a malformed-header, I/O or memory error.
    pub fn open(
        path: impl AsRef<Path>,
        memory_limit_bytes: usize,
    ) -> Result<Self, TrajectoryError> {
        if memory_limit_bytes < BUFFER + size_of::<Self>() {
            return Err(memory(BUFFER + size_of::<Self>(), memory_limit_bytes));
        }
        let file = File::open(path).map_err(|error| io_error(&error))?;
        let file_len = file.metadata().map_err(|error| io_error(&error))?.len();
        let mut input = BufReader::with_capacity(BUFFER, file);
        let pending = next_header(&mut input)?.ok_or_else(invalid)?;
        Ok(Self {
            atoms: pending.atoms,
            input,
            pending: Some(pending),
            next: 0,
            previous_time: None,
            last_step: None,
            last_precision: None,
            limit: memory_limit_bytes,
            file_len,
        })
    }

    /// Simulation step of the most recently emitted frame.
    #[must_use]
    pub const fn last_step(&self) -> Option<i32> {
        self.last_step
    }

    /// Numeric precision of the most recently emitted frame.
    #[must_use]
    pub const fn last_precision(&self) -> Option<TrrPrecision> {
        self.last_precision
    }
}

impl TrajectoryReader for TrrReader {
    fn format(&self) -> &'static str {
        "trr"
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
    fn seek(&mut self, _: usize) -> Result<(), TrajectoryError> {
        Err(TrajectoryError::RandomAccessUnavailable)
    }

    fn read_next(&mut self, frame: &mut Timestep) -> Result<bool, TrajectoryError> {
        self.read_next_bounded(frame, self.limit)
    }

    fn read_next_bounded(
        &mut self,
        frame: &mut Timestep,
        bytes: usize,
    ) -> Result<bool, TrajectoryError> {
        if self.pending.is_none() {
            self.pending = next_header(&mut self.input)?;
        }
        let Some(header) = self.pending.as_ref() else {
            return Ok(false);
        };
        let limit = bytes.min(self.limit);
        let precision = prepare(header, frame, self.atoms, limit)?;
        let payload = header
            .sizes
            .iter()
            .try_fold(0_u64, |sum, size| sum.checked_add(*size as u64))
            .ok_or_else(invalid)?;
        let offset = self
            .input
            .stream_position()
            .map_err(|error| io_error(&error))?;
        if offset
            .checked_add(payload)
            .is_none_or(|end| end > self.file_len)
        {
            return Err(invalid());
        }
        frame.data.clear();
        if header.velocities() == 0 {
            frame.velocities = None;
        }
        if header.forces() == 0 {
            frame.forces = None;
        }
        skip(&mut self.input, header.input_record() + header.energies())?;
        frame.cell = if header.cell() == 0 {
            None
        } else {
            if header.cell() != 9 * precision.bytes() {
                return Err(invalid());
            }
            let mut vectors = [[0.0; 3]; 3];
            for row in &mut vectors {
                for value in row {
                    *value = real(&mut self.input, precision)? * 10.0;
                }
            }
            Some(crate::cell::cell_from_vectors(vectors).ok_or_else(invalid)?)
        };
        skip(
            &mut self.input,
            header.virial() + header.pressure() + header.topology() + header.symmetry(),
        )?;
        fill(
            &mut self.input,
            &mut frame.positions,
            self.atoms,
            precision,
            10.0,
            limit,
        )?;
        if header.velocities() != 0 {
            fill(
                &mut self.input,
                frame.velocities.get_or_insert_with(Vec::new),
                self.atoms,
                precision,
                10.0,
                limit,
            )?;
        }
        if header.forces() != 0 {
            fill(
                &mut self.input,
                frame.forces.get_or_insert_with(Vec::new),
                self.atoms,
                precision,
                0.1,
                limit,
            )?;
        }
        frame.frame = self.next;
        frame.time = Some(header.time);
        frame.dt = self.previous_time.map(|previous| header.time - previous);
        self.previous_time = Some(header.time);
        self.last_step = Some(header.step);
        self.last_precision = Some(precision);
        self.pending = None;
        if self.next == 0 {
            self.pending = next_header(&mut self.input)?;
            frame.dt = self
                .pending
                .as_ref()
                .map(|next| next.time - frame.time.map_or(next.time, |time| time));
        }
        self.next += 1;
        Ok(true)
    }
}

fn prepare(
    header: &Header,
    frame: &Timestep,
    atoms: usize,
    limit: usize,
) -> Result<TrrPrecision, TrajectoryError> {
    if header.atoms != atoms {
        return Err(invalid());
    }
    let precision = header.precision().map_err(|_| invalid())?;
    let atom_bytes = atoms
        .checked_mul(3)
        .and_then(|n| n.checked_mul(precision.bytes()))
        .ok_or_else(invalid)?;
    if header.positions() != atom_bytes {
        return Err(invalid());
    }
    for size in &header.sizes[7..] {
        if *size != 0 && *size != atom_bytes {
            return Err(invalid());
        }
    }
    let vectors = frame
        .positions
        .capacity()
        .max(atoms)
        .checked_add(if header.velocities() != 0 {
            frame
                .velocities
                .as_ref()
                .map_or(atoms, |v| v.capacity().max(atoms))
        } else {
            0
        })
        .and_then(|n| {
            n.checked_add(if header.forces() != 0 {
                frame
                    .forces
                    .as_ref()
                    .map_or(atoms, |v| v.capacity().max(atoms))
            } else {
                0
            })
        })
        .ok_or_else(invalid)?;
    let required = vectors
        .checked_mul(12)
        .and_then(|n| n.checked_add(BUFFER + size_of::<TrrReader>() + size_of::<Timestep>()))
        .ok_or_else(invalid)?;
    if required > limit {
        return Err(memory(required, limit));
    }
    Ok(precision)
}

fn next_header(input: &mut BufReader<File>) -> Result<Option<Header>, TrajectoryError> {
    let mut bytes = [0_u8; 92];
    if input
        .read(&mut bytes[..1])
        .map_err(|error| io_error(&error))?
        == 0
    {
        return Ok(None);
    }
    input
        .read_exact(&mut bytes[1..76])
        .map_err(|error| io_error(&error))?;
    // Header parsing over zero time fields determines precision without allocating.
    let probe = read_header(&mut Xdr::new(&bytes)).map_err(|_| invalid())?;
    let end = 76 + 2 * probe.precision().map_err(|_| invalid())?.bytes();
    input
        .read_exact(&mut bytes[76..end])
        .map_err(|error| io_error(&error))?;
    read_header(&mut Xdr::new(&bytes[..end]))
        .map(Some)
        .map_err(|_| invalid())
}

fn real(input: &mut BufReader<File>, precision: TrrPrecision) -> Result<f64, TrajectoryError> {
    let mut bytes = [0_u8; 8];
    input
        .read_exact(&mut bytes[..precision.bytes()])
        .map_err(|error| io_error(&error))?;
    Xdr::new(&bytes).real(precision).map_err(|_| invalid())
}

fn fill(
    input: &mut BufReader<File>,
    values: &mut Vec<[f32; 3]>,
    atoms: usize,
    precision: TrrPrecision,
    scale: f32,
    limit: usize,
) -> Result<(), TrajectoryError> {
    if values.capacity() < atoms {
        values
            .try_reserve_exact(atoms.saturating_sub(values.len()))
            .map_err(|_| memory(usize::MAX, limit))?;
    }
    values.clear();
    for _ in 0..atoms {
        let mut point = [0.0; 3];
        for value in &mut point {
            *value =
                crate::numeric::f32_from_f64(real(input, precision)?).ok_or_else(invalid)? * scale;
            if !value.is_finite() {
                return Err(invalid());
            }
        }
        values.push(point);
    }
    Ok(())
}

fn skip(input: &mut BufReader<File>, bytes: usize) -> Result<(), TrajectoryError> {
    input
        .seek_relative(i64::try_from(bytes).map_err(|_| invalid())?)
        .map_err(|error| io_error(&error))
}
fn invalid() -> TrajectoryError {
    TrajectoryError::InvalidSource { format: "trr" }
}
fn io_error(error: &std::io::Error) -> TrajectoryError {
    TrajectoryError::SourceIo {
        format: "trr",
        kind: error.kind(),
    }
}
fn memory(required: usize, limit: usize) -> TrajectoryError {
    TrajectoryError::MemoryLimit { required, limit }
}

#[cfg(test)]
#[path = "trr_reader_tests.rs"]
mod tests;

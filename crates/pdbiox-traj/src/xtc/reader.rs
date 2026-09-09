//! Bounded, pull-based XTC file reader.

use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::mem::size_of;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;

use crate::{Timestep, TrajectoryReaderOptions};

use super::{NM_TO_ANGSTROM, XtcError, cell_from_box};

#[path = "reader_support.rs"]
mod support;

#[path = "reader_trait.rs"]
mod reader_trait;

use support::{BoundedBuffer, FrameLayout, guarded_record, memory_error, reserve_exact};

const INPUT_BUFFER_BYTES: usize = 64 * 1024;
const OFFSET_INITIAL_CAPACITY: usize = 1_024;
const FORMAT: &str = "xtc";

/// File-backed XTC source that retains only one decoded frame and bounded scratch.
#[derive(Debug)]
pub struct XtcReader {
    source: molly::XTCReader<BufReader<File>>,
    source_frame: molly::Frame,
    scratch: Vec<u8>,
    offsets: Option<Vec<u64>>,
    file_len: u64,
    next_offset: u64,
    next_frame: usize,
    atom_count: usize,
    frame_count: Option<usize>,
    previous_time: Option<f64>,
    last_step: Option<u32>,
    last_precision: Option<f32>,
    memory_limit_bytes: usize,
    index_unavailable: bool,
}

impl XtcReader {
    /// Opens an XTC file under the default 100 MB reader ceiling.
    ///
    /// # Errors
    ///
    /// Returns an I/O, header, allocation or memory-ceiling error.
    pub fn open(path: &Path) -> Result<Self, XtcError> {
        Self::open_with_memory_limit(path, TrajectoryReaderOptions::DEFAULT_MEMORY_LIMIT_BYTES)
    }

    /// Opens an XTC file under an explicit ceiling of at most 500 MB.
    ///
    /// The ceiling includes the input buffer, decoder storage, lazy offset
    /// index and one caller-owned position buffer.
    ///
    /// # Errors
    ///
    /// Returns an I/O, header, allocation or invalid-memory-limit error.
    pub fn open_with_memory_limit(
        path: &Path,
        memory_limit_bytes: usize,
    ) -> Result<Self, XtcError> {
        if memory_limit_bytes == 0 {
            return Err(XtcError::InvalidMemoryLimit {
                requested: memory_limit_bytes,
            });
        }
        let file = File::open(path)?;
        let file_len = file.metadata()?.len();
        if file_len == 0 {
            return Err(XtcError::InvalidFrame);
        }
        let input = BufReader::with_capacity(INPUT_BUFFER_BYTES, file);
        let mut source = molly::XTCReader::new(input);
        let header = guarded_record(|| source.read_header())?;
        if header.natoms == 0 {
            return Err(XtcError::InvalidFrame);
        }
        source.file.seek(SeekFrom::Start(0))?;
        let scalar_count = header
            .natoms
            .checked_mul(3)
            .ok_or(memory_error(usize::MAX, memory_limit_bytes))?;
        let output_bytes = header
            .natoms
            .checked_mul(size_of::<[f32; 3]>())
            .ok_or(memory_error(usize::MAX, memory_limit_bytes))?;
        let minimum = source
            .file
            .capacity()
            .checked_add(output_bytes.saturating_mul(2))
            .ok_or(memory_error(usize::MAX, memory_limit_bytes))?;
        if minimum > memory_limit_bytes {
            return Err(memory_error(minimum, memory_limit_bytes));
        }
        let mut source_frame = molly::Frame::default();
        source_frame
            .positions
            .try_reserve_exact(scalar_count)
            .map_err(|_| memory_error(usize::MAX, memory_limit_bytes))?;
        let reader = Self {
            source,
            source_frame,
            scratch: Vec::new(),
            offsets: None,
            file_len,
            next_offset: 0,
            next_frame: 0,
            atom_count: header.natoms,
            frame_count: None,
            previous_time: None,
            last_step: None,
            last_precision: None,
            memory_limit_bytes,
            index_unavailable: false,
        };
        let required = reader.total_bytes(header.natoms, 0, 0);
        if required > memory_limit_bytes {
            return Err(memory_error(required, memory_limit_bytes));
        }
        Ok(reader)
    }

    /// Reads one frame into caller-owned storage, retaining its allocation.
    ///
    /// # Errors
    ///
    /// Returns an I/O, malformed-record, allocation or memory-ceiling error.
    pub fn read_next_frame(&mut self, timestep: &mut Timestep) -> Result<bool, XtcError> {
        if self.next_offset == self.file_len {
            self.frame_count = Some(self.next_frame);
            return Ok(false);
        }
        if self.next_offset > self.file_len {
            return Err(XtcError::InvalidFrame);
        }
        let layout = self.read_layout_current(self.next_offset)?;
        if let Err(error) = self.prepare_buffers(timestep, layout.scratch_bytes) {
            self.source.file.seek(SeekFrom::Start(self.next_offset))?;
            return Err(error);
        }
        self.decode(&layout)?;
        let time = f64::from(self.source_frame.time);
        let precision = self.source_frame.precision;
        if !time.is_finite()
            || !precision.is_finite()
            || precision <= 0.0
            || self.source_frame.positions.len() != self.atom_count.saturating_mul(3)
            || self
                .source_frame
                .positions
                .iter()
                .any(|value| !value.is_finite() || !(*value * NM_TO_ANGSTROM).is_finite())
        {
            return Err(XtcError::InvalidFrame);
        }
        let cell = cell_from_box(self.source_frame.boxvec)?;
        let dt = self
            .previous_time
            .map(|previous| time - previous)
            .filter(|value| value.is_finite());
        timestep.positions.clear();
        timestep
            .positions
            .extend(self.source_frame.positions.chunks_exact(3).map(|value| {
                [
                    value[0] * NM_TO_ANGSTROM,
                    value[1] * NM_TO_ANGSTROM,
                    value[2] * NM_TO_ANGSTROM,
                ]
            }));
        timestep.frame = self.next_frame;
        timestep.time = Some(time);
        timestep.dt = dt;
        timestep.velocities = None;
        timestep.forces = None;
        timestep.cell = cell;
        timestep.data.clear();
        self.last_step = Some(self.source_frame.step);
        self.last_precision = Some(precision);
        self.previous_time = Some(time);
        self.next_offset = layout.end;
        self.next_frame = self
            .next_frame
            .checked_add(1)
            .ok_or(XtcError::InvalidFrame)?;
        self.source.step = self.next_frame;
        Ok(true)
    }

    /// Simulation step of the last successfully published frame.
    #[must_use]
    pub const fn last_step(&self) -> Option<u32> {
        self.last_step
    }

    /// Quantization of the last successfully published frame.
    #[must_use]
    pub const fn last_precision(&self) -> Option<f32> {
        self.last_precision
    }

    /// Reader-owned allocation capacity, excluding the caller's timestep.
    #[must_use]
    pub fn workspace_bytes(&self) -> usize {
        self.reader_bytes(self.scratch.capacity(), self.index_capacity())
    }

    fn read_layout_current(&mut self, offset: u64) -> Result<FrameLayout, XtcError> {
        let header = guarded_record(|| self.source.read_header())?;
        if header.natoms != self.atom_count {
            return Err(XtcError::InvalidFrame);
        }
        let positions_offset = offset
            .checked_add(u64::try_from(molly::Header::SIZE).map_err(|_| XtcError::InvalidFrame)?)
            .ok_or(XtcError::InvalidFrame)?;
        let (end, scratch_bytes) = if header.natoms <= 9 {
            let bytes = header
                .natoms
                .checked_mul(size_of::<[f32; 3]>())
                .ok_or(XtcError::InvalidFrame)?;
            let end = positions_offset
                .checked_add(u64::try_from(bytes).map_err(|_| XtcError::InvalidFrame)?)
                .ok_or(XtcError::InvalidFrame)?;
            (end, 0)
        } else {
            self.source.file.seek_relative(32)?;
            let count = molly::reader::read_nbytes(&mut self.source.file, header.magic)
                .map_err(support::record_io)?;
            let padded = count
                .checked_add(molly::padding(count))
                .ok_or(XtcError::InvalidFrame)?;
            let count_bytes = match header.magic {
                molly::Magic::Xtc1995 => 4_u64,
                molly::Magic::Xtc2023 => 8_u64,
            };
            let payload = positions_offset
                .checked_add(32)
                .and_then(|offset| offset.checked_add(count_bytes))
                .ok_or(XtcError::InvalidFrame)?;
            let end = payload
                .checked_add(u64::try_from(padded).map_err(|_| XtcError::InvalidFrame)?)
                .ok_or(XtcError::InvalidFrame)?;
            let rewind = i64::try_from(32 + count_bytes).map_err(|_| XtcError::InvalidFrame)?;
            self.source.file.seek_relative(-rewind)?;
            (end, padded)
        };
        if end > self.file_len || end <= offset {
            return Err(XtcError::InvalidFrame);
        }
        let skip_bytes = end
            .checked_sub(positions_offset)
            .ok_or(XtcError::InvalidFrame)?;
        Ok(FrameLayout {
            header,
            end,
            skip_bytes,
            scratch_bytes,
        })
    }

    fn prepare_buffers(
        &mut self,
        timestep: &mut Timestep,
        scratch_bytes: usize,
    ) -> Result<(), XtcError> {
        let output_capacity = timestep.positions.capacity().max(self.atom_count);
        let scratch_capacity = self.scratch.capacity().max(scratch_bytes);
        let required = self.total_bytes(output_capacity, scratch_capacity, self.index_capacity());
        if required > self.memory_limit_bytes {
            return Err(memory_error(required, self.memory_limit_bytes));
        }
        reserve_exact(
            &mut timestep.positions,
            self.atom_count,
            self.memory_limit_bytes,
        )?;
        // Reserve only. The decoder clears the buffer and sizes it itself, so
        // filling it here would zero the payload twice per frame — on a
        // million-atom trajectory that is megabytes of pointless memset in the
        // frame loop.
        reserve_exact(&mut self.scratch, scratch_bytes, self.memory_limit_bytes)?;
        let actual = self.total_bytes(
            timestep.positions.capacity(),
            self.scratch.capacity(),
            self.index_capacity(),
        );
        if actual > self.memory_limit_bytes {
            return Err(memory_error(actual, self.memory_limit_bytes));
        }
        Ok(())
    }

    fn decode(&mut self, layout: &FrameLayout) -> Result<(), XtcError> {
        self.source_frame.precision = 1_000.0;
        let decoded = catch_unwind(AssertUnwindSafe(|| {
            if layout.header.natoms <= 9 {
                self.source.read_smol_positions(
                    layout.header.natoms,
                    &mut self.source_frame,
                    &molly::selection::AtomSelection::All,
                )
            } else {
                molly::read_positions::<BoundedBuffer<'_>, _>(
                    &mut self.source.file,
                    layout.header.natoms,
                    &mut self.scratch,
                    &mut self.source_frame,
                    &molly::selection::AtomSelection::All,
                    layout.header.magic,
                )
            }
        }))
        .map_err(|_| XtcError::InconsistentRecord)?;
        decoded.map_err(support::record_io)?;
        self.source_frame.step = layout.header.step;
        self.source_frame.time = layout.header.time;
        self.source_frame.boxvec = layout.header.boxvec;
        Ok(())
    }

    fn ensure_index(&mut self) -> Result<(), XtcError> {
        if self.offsets.is_some() {
            return Ok(());
        }
        if self.index_unavailable {
            return Err(memory_error(
                self.memory_limit_bytes.saturating_add(1),
                self.memory_limit_bytes,
            ));
        }
        let saved_offset = self.next_offset;
        let mut offsets = Vec::new();
        let mut cursor = 0_u64;
        let mut max_scratch = self.scratch.capacity();
        self.source.file.seek(SeekFrom::Start(0))?;
        let result = (|| {
            while cursor < self.file_len {
                let layout = self.read_layout_current(cursor)?;
                max_scratch = max_scratch.max(layout.scratch_bytes);
                self.reserve_index_slot(&mut offsets, max_scratch)?;
                offsets.push(cursor);
                let skip = i64::try_from(layout.skip_bytes).map_err(|_| XtcError::InvalidFrame)?;
                self.source.file.seek_relative(skip)?;
                cursor = layout.end;
            }
            reserve_exact(&mut self.scratch, max_scratch, self.memory_limit_bytes)?;
            let required =
                self.total_bytes(self.atom_count, self.scratch.capacity(), offsets.capacity());
            if required > self.memory_limit_bytes {
                return Err(memory_error(required, self.memory_limit_bytes));
            }
            Ok(())
        })();
        self.source.file.seek(SeekFrom::Start(saved_offset))?;
        if let Err(error) = result {
            if matches!(error, XtcError::MemoryLimit { .. }) {
                self.index_unavailable = true;
            }
            return Err(error);
        }
        self.frame_count = Some(offsets.len());
        self.offsets = Some(offsets);
        Ok(())
    }

    fn reserve_index_slot(
        &self,
        offsets: &mut Vec<u64>,
        scratch_capacity: usize,
    ) -> Result<(), XtcError> {
        if offsets.len() < offsets.capacity() {
            let required = self.total_bytes(self.atom_count, scratch_capacity, offsets.capacity());
            return (required <= self.memory_limit_bytes)
                .then_some(())
                .ok_or(memory_error(required, self.memory_limit_bytes));
        }
        let fixed = self.total_bytes(self.atom_count, scratch_capacity, 0);
        let available = self.memory_limit_bytes.saturating_sub(fixed) / size_of::<u64>();
        if offsets.len() >= available {
            return Err(memory_error(
                fixed.saturating_add((offsets.len() + 1).saturating_mul(size_of::<u64>())),
                self.memory_limit_bytes,
            ));
        }
        let growth = offsets.capacity().max(OFFSET_INITIAL_CAPACITY);
        let target = offsets.len().saturating_add(growth).min(available);
        offsets
            .try_reserve_exact(target.saturating_sub(offsets.len()))
            .map_err(|_| memory_error(usize::MAX, self.memory_limit_bytes))?;
        let required = self.total_bytes(self.atom_count, scratch_capacity, offsets.capacity());
        if required > self.memory_limit_bytes {
            return Err(memory_error(required, self.memory_limit_bytes));
        }
        Ok(())
    }

    fn reader_bytes(&self, scratch_capacity: usize, index_capacity: usize) -> usize {
        size_of::<Self>()
            .saturating_add(self.source.file.capacity())
            .saturating_add(
                self.source_frame
                    .positions
                    .capacity()
                    .saturating_mul(size_of::<f32>()),
            )
            .saturating_add(scratch_capacity)
            .saturating_add(index_capacity.saturating_mul(size_of::<u64>()))
    }

    fn total_bytes(
        &self,
        output_capacity: usize,
        scratch_capacity: usize,
        index_capacity: usize,
    ) -> usize {
        self.reader_bytes(scratch_capacity, index_capacity)
            .saturating_add(size_of::<Timestep>())
            .saturating_add(output_capacity.saturating_mul(size_of::<[f32; 3]>()))
    }

    fn index_capacity(&self) -> usize {
        self.offsets.as_ref().map_or(0, Vec::capacity)
    }
}

#[cfg(test)]
#[path = "reader_tests.rs"]
mod tests;

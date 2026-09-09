//! Allocation-bounded adapter around molly's decompression buffer.

use std::io::{self, Read};
use std::panic::{AssertUnwindSafe, catch_unwind};

use crate::{TrajectoryError, XtcError};

pub(super) struct FrameLayout {
    pub(super) header: molly::Header,
    pub(super) end: u64,
    pub(super) skip_bytes: u64,
    pub(super) scratch_bytes: usize,
}

pub(super) struct BoundedBuffer<'a> {
    head: usize,
    scratch: &'a [u8],
}

impl<'s, 'r, R: Read> molly::buffer::Buffered<'s, 'r, R> for BoundedBuffer<'s> {
    /// The decoder clears the scratch buffer immediately before handing it
    /// over, so the preflighted payload size cannot be carried in its length.
    /// Capacity survives that clear, and capacity is what the preflight
    /// actually reserved and charged against the caller's ceiling — so that is
    /// what the payload is checked against. A frame whose declared payload
    /// outgrows the reservation is refused instead of growing the buffer
    /// mid-stream, which is the property the memory contract depends on.
    fn new(scratch: &'s mut Vec<u8>, reader: &'r mut R, magic: molly::Magic) -> io::Result<Self> {
        let count = molly::reader::read_nbytes(reader, magic)?;
        let padded = count
            .checked_add(molly::padding(count))
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "XTC payload overflow"))?;
        if padded > scratch.capacity() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "XTC payload exceeds the bounded preflight",
            ));
        }
        scratch.clear();
        scratch.resize(padded, 0);
        reader.read_exact(scratch)?;
        Ok(Self { head: 0, scratch })
    }

    fn pop(&mut self) -> u8 {
        assert!(!self.scratch.is_empty(), "XTC decoder exceeded its payload");
        let byte = self.scratch[0];
        self.scratch = &self.scratch[1..];
        self.head += 1;
        byte
    }

    fn tell(&self) -> usize {
        self.head
    }

    fn finish(self) -> io::Result<()> {
        Ok(())
    }
}

pub(super) fn reserve_exact<T>(
    values: &mut Vec<T>,
    len: usize,
    limit: usize,
) -> Result<(), XtcError> {
    if len > values.capacity() {
        values
            .try_reserve_exact(len.saturating_sub(values.len()))
            .map_err(|_| memory_error(usize::MAX, limit))?;
    }
    Ok(())
}

pub(super) fn guarded_record<T>(read: impl FnOnce() -> io::Result<T>) -> Result<T, XtcError> {
    catch_unwind(AssertUnwindSafe(read))
        .map_err(|_| XtcError::InconsistentRecord)?
        .map_err(record_io)
}

pub(super) fn record_io(error: io::Error) -> XtcError {
    match error.kind() {
        io::ErrorKind::UnexpectedEof => XtcError::InvalidFrame,
        io::ErrorKind::InvalidData | io::ErrorKind::Other => XtcError::InconsistentRecord,
        _ => XtcError::Io(error),
    }
}

pub(super) const fn memory_error(required: usize, limit: usize) -> XtcError {
    XtcError::MemoryLimit { required, limit }
}

pub(super) fn trajectory_error(error: XtcError) -> TrajectoryError {
    match error {
        XtcError::Io(error) => TrajectoryError::SourceIo {
            format: "xtc",
            kind: error.kind(),
        },
        XtcError::MemoryLimit { required, limit } => {
            TrajectoryError::MemoryLimit { required, limit }
        }
        // The only invalid limit is zero, so the shortfall is the one byte the
        // reader needs before it can hold anything at all.
        XtcError::InvalidMemoryLimit { requested } => TrajectoryError::MemoryLimit {
            required: 1,
            limit: requested,
        },
        XtcError::InvalidFrame | XtcError::InconsistentRecord => {
            TrajectoryError::InvalidSource { format: "xtc" }
        }
    }
}

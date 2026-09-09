//! Window-reusing `MessagePack` cursor with constant payload-skip memory.

use pdbiox_core::{Code, Diagnostic, ExecutionContext, MemoryReservation, SourceBytes};

const MAX_NAME_BYTES: usize = 4_096;

pub(super) struct Cursor<'a, S> {
    source: &'a mut S,
    offset: u64,
    cache_start: u64,
    cache: Vec<u8>,
    window_bytes: usize,
    _reservation: MemoryReservation,
}

impl<'a, S: SourceBytes> Cursor<'a, S> {
    pub(super) fn new(
        source: &'a mut S,
        window_bytes: usize,
        context: &ExecutionContext,
    ) -> Result<Self, Diagnostic> {
        if window_bytes == 0 {
            return Err(error("BinaryCIF index window must not be empty"));
        }
        let reservation = context.try_reserve(window_bytes).map_err(memory_error)?;
        Ok(Self {
            source,
            offset: 0,
            cache_start: 0,
            cache: Vec::with_capacity(window_bytes),
            window_bytes,
            _reservation: reservation,
        })
    }

    pub(super) const fn offset(&self) -> u64 {
        self.offset
    }

    pub(super) fn byte(&mut self) -> Result<u8, Diagnostic> {
        self.ensure_byte()?;
        let local = self.local_offset()?;
        let Some(value) = self.cache.get(local).copied() else {
            return Err(error("truncated BinaryCIF MessagePack value"));
        };
        self.offset = self.checked_advance(1)?;
        Ok(value)
    }

    pub(super) fn unsigned(&mut self, bytes: usize) -> Result<u64, Diagnostic> {
        let mut value = 0_u64;
        for _ in 0..bytes {
            value = (value << 8) | u64::from(self.byte()?);
        }
        Ok(value)
    }

    pub(super) fn text(&mut self, length: u64) -> Result<String, Diagnostic> {
        let length = usize::try_from(length).map_err(|_| error("MessagePack name is too long"))?;
        if length > MAX_NAME_BYTES {
            return Err(error("MessagePack category or item name is too long"));
        }
        let mut bytes = Vec::new();
        bytes
            .try_reserve_exact(length)
            .map_err(|failure| error(failure.to_string()))?;
        for _ in 0..length {
            bytes.push(self.byte()?);
        }
        String::from_utf8(bytes).map_err(|failure| error(failure.to_string()))
    }

    pub(super) fn advance(&mut self, length: u64) -> Result<(), Diagnostic> {
        self.offset = self.checked_advance(length)?;
        Ok(())
    }

    fn ensure_byte(&mut self) -> Result<(), Diagnostic> {
        let cache_length = u64::try_from(self.cache.len())
            .map_err(|_| error("BinaryCIF index cache length is unrepresentable"))?;
        let cache_end = self
            .cache_start
            .checked_add(cache_length)
            .ok_or_else(|| error("BinaryCIF index cache offset overflow"))?;
        if self.offset >= self.cache_start && self.offset < cache_end {
            return Ok(());
        }
        let window = self.source.window(self.offset, self.window_bytes)?;
        self.cache.clear();
        self.cache.extend_from_slice(window.bytes());
        self.cache_start = self.offset;
        Ok(())
    }

    fn local_offset(&self) -> Result<usize, Diagnostic> {
        let delta = self
            .offset
            .checked_sub(self.cache_start)
            .ok_or_else(|| error("BinaryCIF cache precedes its requested offset"))?;
        usize::try_from(delta).map_err(|_| error("BinaryCIF cache offset is unrepresentable"))
    }

    fn checked_advance(&self, length: u64) -> Result<u64, Diagnostic> {
        self.offset
            .checked_add(length)
            .ok_or_else(|| error("BinaryCIF MessagePack offset overflow"))
    }
}

pub(super) fn error(message: impl Into<Box<str>>) -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message(message)
}

fn memory_error(error: impl std::fmt::Display) -> Diagnostic {
    Diagnostic::new(Code::E1902).with_context("reason", error.to_string())
}

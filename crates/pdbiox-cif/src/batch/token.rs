//! Reusable token carry over stable source windows.

use pdbiox_core::{
    Code, Diagnostic, ExecutionContext, MemoryReservation, SourceBytes, StructureBatchError,
};

#[derive(Debug)]
pub(super) struct TokenCursor<S> {
    source: S,
    cursor: u64,
    max_token_bytes: usize,
    cache: Vec<u8>,
    cache_start: u64,
    line_start: bool,
    in_comment: bool,
    _reservation: MemoryReservation,
}

#[derive(Clone, Copy)]
pub(super) struct TokenCheckpoint {
    cursor: u64,
    line_start: bool,
    in_comment: bool,
}

impl<S: SourceBytes> TokenCursor<S> {
    pub(super) fn new(
        source: S,
        max_token_bytes: usize,
        context: &ExecutionContext,
    ) -> Result<Self, StructureBatchError> {
        if max_token_bytes == 0 {
            return Err(StructureBatchError::DemandTooSmall {
                required: 1,
                available: 0,
            });
        }
        let reservation = context.try_reserve(max_token_bytes)?;
        let mut cache = Vec::new();
        cache.try_reserve_exact(max_token_bytes).map_err(|error| {
            StructureBatchError::Diagnostic(
                Diagnostic::new(Code::E1901).with_context("reason", error.to_string()),
            )
        })?;
        Ok(Self {
            source,
            cursor: 0,
            max_token_bytes,
            cache,
            cache_start: 0,
            line_start: true,
            in_comment: false,
            _reservation: reservation,
        })
    }

    pub(super) fn next_token(&mut self) -> Result<Option<&str>, StructureBatchError> {
        if !self.skip_ignored()? {
            return Ok(None);
        }
        let line_start = self.line_start;
        let (local, end) = self.locate_token(line_start)?;
        self.cursor = advance(self.cursor, end)?;
        self.line_start = false;
        let end = local.saturating_add(end);
        let bytes = self.cache.get(local..end).ok_or_else(identity_overflow)?;
        let text = std::str::from_utf8(bytes).map_err(|_| invalid_text())?;
        Ok(Some(trim_token(text, line_start)))
    }

    pub(super) const fn checkpoint(&self) -> TokenCheckpoint {
        TokenCheckpoint {
            cursor: self.cursor,
            line_start: self.line_start,
            in_comment: self.in_comment,
        }
    }

    pub(super) fn restore(&mut self, checkpoint: TokenCheckpoint) {
        self.cursor = checkpoint.cursor;
        self.line_start = checkpoint.line_start;
        self.in_comment = checkpoint.in_comment;
    }

    fn locate_token(&mut self, line_start: bool) -> Result<(usize, usize), StructureBatchError> {
        self.ensure_cursor()?;
        let local = self.local_cursor()?;
        let finished = self.source_finished()?;
        let bytes = self.cache.get(local..).ok_or_else(identity_overflow)?;
        match token_end(bytes, line_start, finished) {
            Ok(end) => Ok((local, end)),
            Err(StructureBatchError::RecordExceedsBudget { .. }) if local != 0 => {
                self.fill(self.cursor)?;
                let finished = self.source_finished()?;
                let end = token_end(&self.cache, line_start, finished)?;
                Ok((0, end))
            }
            Err(error) => Err(error),
        }
    }

    fn skip_ignored(&mut self) -> Result<bool, StructureBatchError> {
        loop {
            self.ensure_cursor()?;
            let local = self.local_cursor()?;
            let bytes = self.cache.get(local..).ok_or_else(identity_overflow)?;
            if bytes.is_empty() {
                return Ok(false);
            }
            let available = bytes.len();
            let mut consumed = 0usize;
            while let Some(byte) = bytes.get(consumed).copied() {
                if self.in_comment {
                    let Some(end) = bytes[consumed..].iter().position(|value| *value == b'\n')
                    else {
                        consumed = bytes.len();
                        break;
                    };
                    consumed = consumed.saturating_add(end + 1);
                    self.in_comment = false;
                    self.line_start = true;
                    continue;
                }
                if byte.is_ascii_whitespace() {
                    self.line_start = byte == b'\n';
                    consumed += 1;
                    continue;
                }
                if byte == b'#' {
                    self.in_comment = true;
                    continue;
                }
                break;
            }
            self.cursor = advance(self.cursor, consumed)?;
            if consumed < available {
                return Ok(true);
            }
        }
    }

    fn ensure_cursor(&mut self) -> Result<(), StructureBatchError> {
        let cache_length = u64::try_from(self.cache.len()).map_err(|_| identity_overflow())?;
        let cache_end = self
            .cache_start
            .checked_add(cache_length)
            .ok_or_else(identity_overflow)?;
        if self.cache.is_empty() || self.cursor < self.cache_start || self.cursor >= cache_end {
            self.fill(self.cursor)?;
        }
        Ok(())
    }

    fn fill(&mut self, start: u64) -> Result<(), StructureBatchError> {
        let window = self.source.window(start, self.max_token_bytes)?;
        self.cache.clear();
        self.cache.extend_from_slice(window.bytes());
        self.cache_start = start;
        Ok(())
    }

    fn local_cursor(&self) -> Result<usize, StructureBatchError> {
        let local = self
            .cursor
            .checked_sub(self.cache_start)
            .ok_or_else(identity_overflow)?;
        usize::try_from(local).map_err(|_| identity_overflow())
    }

    fn source_finished(&self) -> Result<bool, StructureBatchError> {
        let Some(length) = self.source.len_hint() else {
            return Ok(self.cache.len() < self.max_token_bytes);
        };
        let cache_length = u64::try_from(self.cache.len()).map_err(|_| identity_overflow())?;
        let cache_end = self
            .cache_start
            .checked_add(cache_length)
            .ok_or_else(identity_overflow)?;
        Ok(cache_end >= length)
    }
}

fn token_end(
    bytes: &[u8],
    line_start: bool,
    source_finished: bool,
) -> Result<usize, StructureBatchError> {
    let Some(first) = bytes.first().copied() else {
        return Ok(0);
    };
    if line_start && first == b';' {
        return terminated(bytes, source_finished);
    }
    if matches!(first, b'\'' | b'"') {
        for position in 1..bytes.len() {
            if bytes[position] == first
                && bytes.get(position + 1).is_none_or(u8::is_ascii_whitespace)
            {
                return Ok(position + 1);
            }
        }
        return incomplete(bytes.len());
    }
    if let Some(end) = bytes.iter().position(u8::is_ascii_whitespace) {
        return Ok(end);
    }
    if source_finished {
        return Ok(bytes.len());
    }
    incomplete(bytes.len())
}

fn terminated(bytes: &[u8], source_finished: bool) -> Result<usize, StructureBatchError> {
    if let Some(position) = bytes.windows(2).skip(1).position(|window| window == b"\n;") {
        return Ok(position.saturating_add(3));
    }
    if source_finished {
        return Err(StructureBatchError::Diagnostic(
            Diagnostic::new(Code::E1102).with_message("multiline CIF value was not terminated"),
        ));
    }
    incomplete(bytes.len())
}

fn incomplete(available: usize) -> Result<usize, StructureBatchError> {
    Err(StructureBatchError::RecordExceedsBudget {
        required: available.saturating_add(1),
        available,
    })
}

fn trim_token(text: &str, line_start: bool) -> &str {
    let bytes = text.as_bytes();
    if line_start && bytes.first() == Some(&b';') {
        return match text.get(1..text.len().saturating_sub(2)) {
            Some(value) => value,
            None => "",
        };
    }
    if matches!(bytes.first(), Some(b'\'' | b'"')) {
        return match text.get(1..text.len().saturating_sub(1)) {
            Some(value) => value,
            None => "",
        };
    }
    text
}

fn advance(start: u64, count: usize) -> Result<u64, StructureBatchError> {
    let count = u64::try_from(count).map_err(|_| identity_overflow())?;
    start.checked_add(count).ok_or_else(identity_overflow)
}

fn identity_overflow() -> StructureBatchError {
    StructureBatchError::Diagnostic(
        Diagnostic::new(Code::E1903).with_message("CIF byte offset exceeds u64"),
    )
}

fn invalid_text() -> StructureBatchError {
    StructureBatchError::Diagnostic(
        Diagnostic::new(Code::E1101).with_message("CIF input is not valid UTF-8"),
    )
}

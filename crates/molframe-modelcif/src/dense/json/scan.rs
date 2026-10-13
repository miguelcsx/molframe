//! Allocation-free cursor over JSON text, specialised for numeric arrays.
//!
//! The importer needs to read a megabyte-scale matrix twice (probe, then fill)
//! without building a tree, so values are parsed straight from the borrowed
//! bytes and unrelated members are skipped structurally.

use crate::ConfidenceError;

/// Cursor over a JSON byte slice.
#[derive(Clone, Copy)]
pub(super) struct Scanner<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Scanner<'a> {
    pub(super) const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, pos: 0 }
    }

    pub(super) const fn at(bytes: &'a [u8], pos: usize) -> Self {
        Self { bytes, pos }
    }

    pub(super) const fn offset(&self) -> usize {
        self.pos
    }

    pub(super) const fn fail(&self, expected: &'static str) -> ConfidenceError {
        ConfidenceError::Json {
            offset: self.pos,
            expected,
        }
    }

    pub(super) fn skip_ws(&mut self) {
        while matches!(self.bytes.get(self.pos), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.pos += 1;
        }
    }

    pub(super) fn peek(&mut self) -> Option<u8> {
        self.skip_ws();
        self.bytes.get(self.pos).copied()
    }

    pub(super) fn expect(
        &mut self,
        byte: u8,
        expected: &'static str,
    ) -> Result<(), ConfidenceError> {
        if self.peek() == Some(byte) {
            self.pos += 1;
            Ok(())
        } else {
            Err(self.fail(expected))
        }
    }

    /// Consumes `byte` when it is next.
    pub(super) fn eat(&mut self, byte: u8) -> bool {
        if self.peek() == Some(byte) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    /// Reads a string without escapes and returns its content.
    pub(super) fn string(&mut self) -> Result<&'a [u8], ConfidenceError> {
        self.require(b'"', "a string")?;
        let start = self.pos;
        loop {
            match self.bytes.get(self.pos) {
                Some(b'"') => {
                    let content = &self.bytes[start..self.pos];
                    self.pos += 1;
                    return Ok(content);
                }
                Some(b'\\') => self.pos += 2,
                Some(_) => self.pos += 1,
                None => return Err(self.fail("a closing quote")),
            }
        }
    }

    /// Reads one finite number.
    pub(super) fn number(&mut self) -> Result<f32, ConfidenceError> {
        self.skip_ws();
        let start = self.pos;
        while matches!(
            self.bytes.get(self.pos),
            Some(b'0'..=b'9' | b'-' | b'+' | b'.' | b'e' | b'E')
        ) {
            self.pos += 1;
        }
        let text =
            std::str::from_utf8(&self.bytes[start..self.pos]).map_err(|_| self.fail("a number"))?;
        match text.parse::<f32>() {
            Ok(value) if value.is_finite() => Ok(value),
            _ => Err(Self::at(self.bytes, start).fail("a finite number")),
        }
    }

    /// Skips one value of any kind.
    pub(super) fn skip_value(&mut self) -> Result<(), ConfidenceError> {
        match self.peek() {
            Some(b'"') => self.string().map(drop),
            Some(b'[' | b'{') => self.skip_container(),
            Some(_) => {
                while !matches!(
                    self.bytes.get(self.pos),
                    None | Some(b',' | b']' | b'}' | b' ' | b'\t' | b'\n' | b'\r')
                ) {
                    self.pos += 1;
                }
                Ok(())
            }
            None => Err(self.fail("a value")),
        }
    }

    fn skip_container(&mut self) -> Result<(), ConfidenceError> {
        let mut depth = 0_usize;
        loop {
            match self.bytes.get(self.pos) {
                Some(b'[' | b'{') => depth += 1,
                Some(b']' | b'}') => depth -= 1,
                Some(b'"') => {
                    self.string()?;
                    if depth == 0 {
                        return Ok(());
                    }
                    continue;
                }
                Some(_) => {}
                None => return Err(self.fail("a closing bracket")),
            }
            self.pos += 1;
            if depth == 0 {
                return Ok(());
            }
        }
    }
}

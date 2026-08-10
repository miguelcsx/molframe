//! Walking a file line by line, tracking where each began.
//!
//! The position is a by-product of the walk rather than something reconstructed
//! afterwards, so attaching an exact location to a finding costs nothing.

use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::span::Position;

/// A line of the file, with where it began.
pub(super) struct Line<'a> {
    pub(super) text: &'a str,
    pub(super) at: Position,
}

/// Walks the lines of a file, tracking where each began.
pub(super) struct Lines<'a> {
    remaining: &'a str,
    at: Position,
}

impl<'a> Lines<'a> {
    pub(super) fn new(text: &'a str) -> Self {
        Self {
            remaining: text,
            at: Position::START,
        }
    }
}

impl<'a> Iterator for Lines<'a> {
    type Item = Result<Line<'a>, Diagnostic>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining.is_empty() {
            return None;
        }
        let at = self.at;
        let (line, rest) = match self.remaining.find('\n') {
            Some(end) => (self.remaining.get(..end)?, self.remaining.get(end + 1..)?),
            None => (self.remaining, ""),
        };
        self.remaining = rest;
        let advance = u32::try_from(line.len())
            .ok()
            .and_then(|length| length.checked_add(1));
        let Some((byte_offset, line_number)) = advance
            .and_then(|advance| at.byte_offset.checked_add(advance))
            .zip(at.line.checked_add(1))
        else {
            self.remaining = "";
            return Some(Err(Diagnostic::new(Code::E1202)
                .with_message("PDB source position exceeds the supported range")));
        };
        self.at = Position::new(byte_offset, line_number, 1);
        Some(Ok(Line {
            text: line.trim_end_matches('\r'),
            at,
        }))
    }
}

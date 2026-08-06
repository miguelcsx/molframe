//! Walking a file line by line, tracking where each began.
//!
//! The position is a by-product of the walk rather than something reconstructed
//! afterwards, so attaching an exact location to a finding costs nothing.

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
    type Item = Line<'a>;

    fn next(&mut self) -> Option<Line<'a>> {
        if self.remaining.is_empty() {
            return None;
        }
        let at = self.at;
        let (line, rest) = match self.remaining.find('\n') {
            Some(end) => (self.remaining.get(..end)?, self.remaining.get(end + 1..)?),
            None => (self.remaining, ""),
        };
        self.remaining = rest;
        self.at = Position::new(at.byte_offset + line.len() as u32 + 1, at.line + 1, 1);
        Some(Line {
            text: line.trim_end_matches('\r'),
            at,
        })
    }
}

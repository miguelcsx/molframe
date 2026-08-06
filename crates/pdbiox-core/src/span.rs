//! Locations within a source file.
//!
//! Every parser carries a cursor that already knows its byte offset, line and
//! column, so attaching an exact location to a diagnostic costs nothing beyond
//! copying three integers. That is the whole reason a hand-written lexer is
//! worth the effort over a combinator stack: the position is a by-product of
//! scanning rather than something reconstructed afterwards.

use std::fmt;

/// A point in a source file.
///
/// `line` and `column` are 1-based, as an editor reports them. `byte_offset` is
/// 0-based, as a slice index requires.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct Position {
    /// Offset from the start of the input, in bytes.
    pub byte_offset: u32,
    /// 1-based line number.
    pub line: u32,
    /// 1-based column number, counted in bytes.
    pub column: u32,
}

impl Position {
    /// The start of a file.
    pub const START: Self = Self {
        byte_offset: 0,
        line: 1,
        column: 1,
    };

    /// Creates a position.
    ///
    /// # Examples
    ///
    /// ```
    /// use pdbiox_core::Position;
    ///
    /// let at = Position::new(42, 3, 7);
    /// assert_eq!(at.line, 3);
    /// ```
    #[must_use]
    pub const fn new(byte_offset: u32, line: u32, column: u32) -> Self {
        Self {
            byte_offset,
            line,
            column,
        }
    }

    /// Advances past one byte, tracking the line and column.
    ///
    /// A newline moves to column 1 of the next line; everything else advances
    /// the column. Carriage returns advance the column like any other byte, so a
    /// CRLF file reports the column of the visible character.
    #[must_use]
    pub const fn advance(self, byte: u8) -> Self {
        let byte_offset = self.byte_offset.saturating_add(1);

        if byte == b'\n' {
            Self {
                byte_offset,
                line: self.line.saturating_add(1),
                column: 1,
            }
        } else {
            Self {
                byte_offset,
                line: self.line,
                column: self.column.saturating_add(1),
            }
        }
    }
}

impl fmt::Display for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.line, self.column)
    }
}

/// A half-open byte range within a source file, with the start's line and column.
///
/// Spans are retained on every category and value a document holds, which is
/// what lets a write that preserves the original file reproduce untouched values
/// byte for byte instead of re-rendering them.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct ByteSpan {
    /// Where the span begins.
    pub start: Position,
    /// Offset one past the last byte of the span.
    pub end: u32,
}

impl ByteSpan {
    /// Creates a span running from `start` up to but excluding `end`.
    ///
    /// # Examples
    ///
    /// ```
    /// use pdbiox_core::{ByteSpan, Position};
    ///
    /// let span = ByteSpan::new(Position::new(10, 2, 1), 14);
    /// assert_eq!(span.len(), 4);
    /// ```
    #[must_use]
    pub const fn new(start: Position, end: u32) -> Self {
        Self { start, end }
    }

    /// Creates an empty span at a point.
    #[must_use]
    pub const fn empty(at: Position) -> Self {
        Self {
            start: at,
            end: at.byte_offset,
        }
    }

    /// Returns the length of the span in bytes.
    #[must_use]
    pub const fn len(self) -> u32 {
        self.end.saturating_sub(self.start.byte_offset)
    }

    /// Returns true when the span covers no bytes.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.len() == 0
    }

    /// Returns the bytes this span covers, or `None` if it runs past the end of
    /// `source`.
    ///
    /// Returning `None` rather than clamping matters: a span that does not fit
    /// its buffer means the span and the buffer came from different reads, and
    /// silently truncating would hide that.
    #[must_use]
    pub fn slice(self, source: &[u8]) -> Option<&[u8]> {
        let start = usize::try_from(self.start.byte_offset).ok()?;
        let end = usize::try_from(self.end).ok()?;

        source.get(start..end)
    }
}

impl fmt::Display for ByteSpan {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.start)
    }
}

#[cfg(test)]
#[path = "span_tests.rs"]
mod tests;

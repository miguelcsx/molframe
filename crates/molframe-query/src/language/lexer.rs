//! Bounded tokenisation for the selection language.

use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::span::{ByteSpan, Position};

/// A lexical token and the exact source range it was read from.
#[derive(Clone, PartialEq, Eq, Debug)]
pub(super) struct Token {
    pub(super) kind: TokenKind,
    pub(super) span: ByteSpan,
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub(super) enum TokenKind {
    Value(Box<str>),
    Operator(Box<str>),
    LeftParen,
    RightParen,
}

/// Tokenizes selection-language source into owned tokens.
///
/// Runtime is `O(N)` in source bytes. Ordinary and unescaped quoted tokens are
/// copied exactly once into their final `Box<str>`. Escaped quoted values use a
/// single exactly-sized temporary `String`. Every token carries its byte range
/// together with the line and column of its first byte; the line cursor
/// advances with the scan, so locating a token costs nothing extra.
pub(super) fn lex(source: &str) -> Result<Vec<Token>, Diagnostic> {
    let bytes = source.as_bytes();
    let mut tokens = Vec::new();
    let mut cursor = Cursor::default();
    let mut position = 0usize;

    while position < bytes.len() {
        if bytes[position].is_ascii_whitespace() {
            cursor.consume(bytes, position, position + 1);
            position += 1;
            continue;
        }

        let start = position;

        let kind = match bytes[position] {
            b'(' => {
                position += 1;
                TokenKind::LeftParen
            }
            b')' => {
                position += 1;
                TokenKind::RightParen
            }
            b'<' | b'>' | b'=' | b'!' => {
                position += 1;

                if position < bytes.len() && bytes[position] == b'=' {
                    position += 1;
                }

                TokenKind::Operator(source[start..position].into())
            }
            b'\'' | b'"' => {
                let (value, next) = quoted_value(source, start, cursor)?;
                position = next;
                TokenKind::Value(value)
            }
            _ => {
                while position < bytes.len()
                    && !bytes[position].is_ascii_whitespace()
                    && !matches!(bytes[position], b'(' | b')' | b'<' | b'>' | b'=' | b'!')
                {
                    position += 1;
                }

                let value = &source[start..position];

                if value.is_empty() {
                    return Err(syntax(
                        cursor.span(start, start + 1),
                        "selection contains an unreadable token",
                    ));
                }

                TokenKind::Value(value.into())
            }
        };

        tokens.push(Token {
            kind,
            span: cursor.span(start, position),
        });
        cursor.consume(bytes, start, position);
    }

    Ok(tokens)
}

/// The position of the scan within the source, tracked incrementally.
#[derive(Clone, Copy, Default)]
pub(super) struct Cursor {
    /// Zero-based line of the next unread byte.
    line: u64,
    /// Byte offset at which that line begins.
    line_start: usize,
}

impl Cursor {
    /// A span over `start..end`, which must not begin before this cursor.
    pub(super) fn span(self, start: usize, end: usize) -> ByteSpan {
        let column = start.saturating_sub(self.line_start);
        ByteSpan::new(
            Position::new(
                widen(start),
                self.line.saturating_add(1),
                widen(column).saturating_add(1),
            ),
            widen(end),
        )
    }

    /// Advances past `bytes[start..end]`, counting any line breaks inside.
    fn consume(&mut self, bytes: &[u8], start: usize, end: usize) {
        let Some(slice) = bytes.get(start..end) else {
            return;
        };
        for (index, byte) in slice.iter().enumerate() {
            if *byte == b'\n' {
                self.line = self.line.saturating_add(1);
                self.line_start = start + index + 1;
            }
        }
    }
}

/// An empty span just past the last byte, where "expected more" is reported.
pub(super) fn end_of(source: &str) -> ByteSpan {
    let mut cursor = Cursor::default();
    cursor.consume(source.as_bytes(), 0, source.len());
    cursor.span(source.len(), source.len())
}

fn widen(value: usize) -> u64 {
    // Lossless: no supported target has a `usize` wider than 64 bits.
    value as u64
}

/// Parses one quoted token beginning at `start`.
///
/// Backslash escapes retain the original lexer semantics: the slash is removed
/// and the immediately following Unicode scalar is copied literally.
fn quoted_value(
    source: &str,
    start: usize,
    cursor: Cursor,
) -> Result<(Box<str>, usize), Diagnostic> {
    let bytes = source.as_bytes();
    let unclosed = || {
        syntax(
            cursor.span(start, source.len()),
            "quoted selection value is not closed",
        )
    };
    let Some(&quote) = bytes.get(start) else {
        return Err(unclosed());
    };

    let Some(content_start) = start.checked_add(1) else {
        return Err(unclosed());
    };
    let mut position = content_start;
    let mut escape_count = 0usize;

    while position < bytes.len() {
        if bytes[position] == quote {
            let value = if escape_count == 0 {
                source[content_start..position].into()
            } else {
                decode_quoted(source, content_start, position, escape_count).ok_or_else(|| {
                    syntax(
                        cursor.span(start, position + 1),
                        "quoted selection bounds are inconsistent",
                    )
                })?
            };

            return Ok((value, position + 1));
        }

        if bytes[position] == b'\\' && position + 1 < bytes.len() {
            escape_count += 1;
            position += 1;
        }

        let Some(character) = source.get(position..).and_then(|tail| tail.chars().next()) else {
            break;
        };

        position += character.len_utf8();
    }

    Err(unclosed())
}

/// Removes lexer escape markers from an already bounded quoted substring.
///
/// Capacity is computed exactly because every escape removes one ASCII byte.
/// Runtime is `O(N)` and requires one allocation.
fn decode_quoted(source: &str, start: usize, end: usize, escape_count: usize) -> Option<Box<str>> {
    let capacity = end
        .checked_sub(start)
        .and_then(|length| length.checked_sub(escape_count))?;
    let mut value = String::with_capacity(capacity);
    let bytes = source.as_bytes();
    let mut position = start;

    while position < end {
        if bytes.get(position) == Some(&b'\\') && position + 1 < end {
            position += 1;
        }

        let Some(character) = source
            .get(position..end)
            .and_then(|tail| tail.chars().next())
        else {
            break;
        };

        value.push(character);
        position += character.len_utf8();
    }

    Some(value.into_boxed_str())
}

/// Creates a lexical syntax diagnostic located at `span`.
///
/// Runtime and auxiliary space are `O(1)` apart from diagnostic-owned context.
pub(super) fn syntax(span: ByteSpan, message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E4001).with_message(message).at(span)
}

#[cfg(test)]
#[path = "lexer_tests.rs"]
mod tests;

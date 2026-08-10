//! Bounded tokenisation for the selection language.

use pdbiox_core::diagnostic::{Code, Diagnostic};

/// A lexical token and its byte offset.
#[derive(Clone, PartialEq, Eq, Debug)]
pub(super) struct Token {
    pub(super) kind: TokenKind,
    pub(super) offset: usize,
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
/// single exactly-sized temporary `String`.
pub(super) fn lex(source: &str) -> Result<Vec<Token>, Diagnostic> {
    let bytes = source.as_bytes();
    let mut tokens = Vec::new();
    let mut position = 0usize;

    while position < bytes.len() {
        if bytes[position].is_ascii_whitespace() {
            position += 1;
            continue;
        }

        let start = position;

        match bytes[position] {
            b'(' => {
                tokens.push(Token {
                    kind: TokenKind::LeftParen,
                    offset: start,
                });
                position += 1;
            }
            b')' => {
                tokens.push(Token {
                    kind: TokenKind::RightParen,
                    offset: start,
                });
                position += 1;
            }
            b'<' | b'>' | b'=' | b'!' => {
                position += 1;

                if position < bytes.len() && bytes[position] == b'=' {
                    position += 1;
                }

                tokens.push(Token {
                    kind: TokenKind::Operator(source[start..position].into()),
                    offset: start,
                });
            }
            b'\'' | b'"' => {
                let (value, next) = quoted_value(source, start)?;

                tokens.push(Token {
                    kind: TokenKind::Value(value),
                    offset: start,
                });

                position = next;
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
                    return Err(syntax(start, "selection contains an unreadable token"));
                }

                tokens.push(Token {
                    kind: TokenKind::Value(value.into()),
                    offset: start,
                });
            }
        }
    }

    Ok(tokens)
}

/// Parses one quoted token beginning at `start`.
///
/// Backslash escapes retain the original lexer semantics: the slash is removed
/// and the immediately following Unicode scalar is copied literally.
fn quoted_value(source: &str, start: usize) -> Result<(Box<str>, usize), Diagnostic> {
    let bytes = source.as_bytes();
    let Some(&quote) = bytes.get(start) else {
        return Err(syntax(start, "quoted selection value is not closed"));
    };

    let Some(content_start) = start.checked_add(1) else {
        return Err(syntax(
            start,
            "quoted selection offset exceeds the addressable range",
        ));
    };
    let mut position = content_start;
    let mut escape_count = 0usize;

    while position < bytes.len() {
        if bytes[position] == quote {
            let value = if escape_count == 0 {
                source[content_start..position].into()
            } else {
                decode_quoted(source, content_start, position, escape_count)?
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

    Err(syntax(start, "quoted selection value is not closed"))
}

/// Removes lexer escape markers from an already bounded quoted substring.
///
/// Capacity is computed exactly because every escape removes one ASCII byte.
/// Runtime is `O(N)` and requires one allocation.
fn decode_quoted(
    source: &str,
    start: usize,
    end: usize,
    escape_count: usize,
) -> Result<Box<str>, Diagnostic> {
    let Some(capacity) = end
        .checked_sub(start)
        .and_then(|length| length.checked_sub(escape_count))
    else {
        return Err(syntax(start, "quoted selection bounds are inconsistent"));
    };
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

    Ok(value.into_boxed_str())
}

/// Creates a lexical syntax diagnostic at byte `offset`.
///
/// Runtime and auxiliary space are `O(1)` apart from diagnostic-owned context.
fn syntax(offset: usize, message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E4001)
        .with_message(message)
        .with_context("offset", offset.to_string())
}

#[cfg(test)]
#[path = "lexer_tests.rs"]
mod tests;

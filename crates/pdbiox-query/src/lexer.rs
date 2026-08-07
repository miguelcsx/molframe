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
                let quote = bytes[position];
                position += 1;
                let mut value = String::new();
                let mut closed = false;
                while position < bytes.len() {
                    if bytes[position] == quote {
                        position += 1;
                        closed = true;
                        break;
                    }
                    if bytes[position] == b'\\' && position + 1 < bytes.len() {
                        position += 1;
                    }
                    let Some(character) = source[position..].chars().next() else {
                        break;
                    };
                    value.push(character);
                    position += character.len_utf8();
                }
                if !closed {
                    return Err(syntax(start, "quoted selection value is not closed"));
                }
                tokens.push(Token {
                    kind: TokenKind::Value(value.into()),
                    offset: start,
                });
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

fn syntax(offset: usize, message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E4001)
        .with_message(message)
        .with_context("offset", offset.to_string())
}

#[cfg(test)]
#[path = "lexer_tests.rs"]
mod tests;
